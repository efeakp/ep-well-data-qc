"""End-to-end audit: discover -> read -> reconcile -> standardise -> condition -> predict -> report."""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import lasio
import numpy as np
import pandas as pd

from . import plots
from .checkshots import qc_checkshot
from .condition import condition_well, depth_match, splice
from .io import (
    read_checkshot,
    read_deviation,
    read_grid_extent,
    read_las,
    read_well_head,
    read_well_tops,
    scan_segy,
    sniff,
)
from .io.segy import same_grid, xy_to_grid
from .model import Issue, Lineage, Well
from .report import curve_scores, issues_frame, write_report
from .rockphysics import dt_to_vp, greenberg_castagna_vs, predict_curve, rpt_distance, rpt_lines, vp_to_dt
from .standardise import (
    assign_master_uwis,
    canonical_name,
    find_close_wells,
    find_duplicate_wells,
    match_by_location,
    match_name,
    standardise_curves,
    to_depth_unit,
)


class Audit:
    """Holds the state of one audit run. Use :func:`run_audit` for the one-call version."""

    def __init__(self, projects: dict[str, Path], depth_unit: str = "m", uwi_prefix: str = "W"):
        self.projects = {k: Path(v) for k, v in projects.items()}
        self.depth_unit = depth_unit
        self.uwi_prefix = uwi_prefix
        self.lineage = Lineage()
        self.issues: list[Issue] = []
        self.inventory: list[dict] = []
        self.wells: list[Well] = []
        self.segy: list[tuple[str, object]] = []
        self.surfaces: list[dict] = []
        self.raw_logs: dict[str, pd.DataFrame] = {}
        self.clean_logs: dict[str, pd.DataFrame] = {}
        self.checkshot_qc: dict[str, pd.DataFrame] = {}
        self.predictions: dict[str, dict[str, pd.Series]] = {}
        self.validation = pd.DataFrame()
        self.blind_predictions: dict[str, pd.Series] = {}
        self.master = pd.DataFrame()
        self.duplicates = pd.DataFrame()
        self.figures: list[str] = []

    # ------------------------------------------------------------------ discovery & reading
    def discover(self) -> None:
        for proj, root in self.projects.items():
            for p in sorted(root.rglob("*")):
                if p.is_file():
                    kind = sniff(p)
                    self.inventory.append({"project": proj, "file": str(p.relative_to(root)), "type": kind,
                                           "size_kb": round(p.stat().st_size / 1024, 1), "path": str(p)})
                    if kind == "unknown":
                        self.issues.append(Issue("minor", "inventory", "File type not recognised",
                                                 source=f"{proj}/{p.name}"))
        self.lineage.add("discover", "inventory files", n_affected=len(self.inventory))

    def _files(self, proj: str, kind: str) -> list[Path]:
        return [Path(r["path"]) for r in self.inventory if r["project"] == proj and r["type"] == kind]

    def read(self) -> None:
        for proj in self.projects:
            heads = []
            for f in self._files(proj, "well_head"):
                h, iss = read_well_head(f)
                heads.append(h)
                self._add(iss, proj)
            heads = pd.concat(heads, ignore_index=True) if heads else pd.DataFrame(
                columns=["name", "uwi", "x", "y", "kb", "td", "xy_unit", "depth_unit"])
            if heads.empty:
                self.issues.append(Issue("major", "inventory", "No well header file in project", source=proj))
            pw: dict[str, Well] = {}
            for _, r in heads.iterrows():
                w = Well(name=str(r["name"]), uwi=r.get("uwi"), x=r["x"], y=r["y"], kb=r["kb"], td=r["td"],
                         depth_unit=r.get("depth_unit") or None, project=proj, sources=[str(r["source"])])
                pw[w.name] = w

            def resolve(name: str | None, uwi: str | None, xy: tuple | None, source: str, kind: str,
                        pw=pw, heads=heads, proj=proj) -> Well:
                names = list(pw)
                m, how = match_name(name, names) if name else (None, "none")
                if m is None and uwi:
                    m = next((n for n, w in pw.items() if w.uwi and str(w.uwi).upper() == str(uwi).upper()), None)
                    how = "uwi" if m else how
                if m is None and xy is not None and not heads.empty:
                    m = match_by_location(xy[0], xy[1], heads)
                    how = "location" if m else how
                if m is None:
                    w = Well(name=name or Path(source).stem, uwi=uwi, project=proj)
                    pw[w.name] = w
                    self.issues.append(Issue("major", "naming", f"{kind} well '{name}' not in well header; "
                                                                "created a new well record", w.name, source=source))
                    return w
                w = pw[m]
                if how in ("fuzzy", "location"):
                    self.issues.append(Issue("major", "naming",
                                             f"{kind} names the well '{name}' but it was matched to '{m}' by "
                                             f"{'similar name' if how == 'fuzzy' else 'wellhead location'}",
                                             m, source=source))
                if name and name != m:
                    w.aliases.add(name)
                return w

            runs: dict[str, list[tuple[str, pd.DataFrame]]] = {}
            for f in self._files(proj, "las"):
                hdr, df, units, iss = read_las(f)
                w = resolve(hdr["well"], hdr["uwi"], None, f.name, "LAS")
                self._add(iss, proj, w.name)
                if hdr["uwi"] and w.uwi and str(hdr["uwi"]) != str(w.uwi):
                    self.issues.append(Issue("major", "naming", f"LAS UWI '{hdr['uwi']}' differs from header "
                                                                f"UWI '{w.uwi}'", w.name, source=f.name))
                runs.setdefault(w.name, []).append((f.name, df))
                w.units = {**w.units, **units}
                lu = hdr["depth_unit"]
                if w.depth_unit and lu and w.depth_unit != lu:
                    self.issues.append(Issue("critical", "units", f"LAS depth unit '{lu}' differs from header "
                                                                  f"'{w.depth_unit}'", w.name, source=f.name))
                w.depth_unit = w.depth_unit or lu
                w.sources.append(f.name)
            for wname, rs in runs.items():
                pw[wname].logs = self._merge_runs(pw[wname], rs)
            for f in self._files(proj, "deviation"):
                meta, traj, iss = read_deviation(f)
                w = resolve(meta.get("name"), None, (meta.get("x"), meta.get("y")), f.name, "Deviation file")
                for i in iss:
                    i.well = w.name
                self._add(iss, proj, w.name)
                w.trajectory = traj
                w.sources.append(f.name)
                self._check_consistency(w, meta, f.name)
            for f in self._files(proj, "well_tops"):
                tops, iss = read_well_tops(f)
                self._add(iss, proj)
                for wn, g in tops.groupby("well"):
                    w = resolve(wn, None, None, f.name, "Well tops")
                    w.tops = pd.concat([w.tops, g]) if w.tops is not None else g
                    w.sources.append(f.name)
            for f in self._files(proj, "checkshot"):
                meta, cs, iss = read_checkshot(f)
                name = meta["well"] or self._name_from_filename(f.stem, list(pw))
                w = resolve(name, None, None, f.name, "Checkshot")
                for i in iss:
                    i.well = i.well or w.name
                self._add(iss, proj, w.name)
                w.checkshot = cs
                w.sources.append(f.name)
            for f in self._files(proj, "segy"):
                info = scan_segy(f)
                self._add(info.issues, proj)
                self.segy.append((proj, info))
            for f in self._files(proj, "surface"):
                g = read_grid_extent(f)
                if g:
                    g["project"] = proj
                    self.surfaces.append(g)
            self.wells.extend(pw.values())
        self.lineage.add("read", "load and reconcile files", n_affected=len(self.wells))

    def _merge_runs(self, w: Well, runs: list[tuple[str, pd.DataFrame]]) -> pd.DataFrame:
        """Depth-match overlapping logging runs on GR, then splice them (shallow run wins in overlaps)."""
        if len(runs) == 1:
            return runs[0][1]
        runs = sorted(runs, key=lambda r: r[1].index.min())
        merged_name, merged = runs[0]
        out_runs = [merged]
        unit = w.depth_unit or "m"
        max_shift = 5.0 if unit == "m" else 15.0
        gr = next((c for c in ("GR", "GRC", "SGR") if c in merged.columns), None)
        for name, df in runs[1:]:
            shift, corr = 0.0, np.nan
            if gr and gr in df.columns:
                try:
                    shift, corr = depth_match(merged[gr], df[gr], max_shift)
                except ValueError:
                    pass
            step = float(np.median(np.diff(df.index)))
            if abs(shift) >= step and corr > 0.6:
                df = df.copy()
                df.index = np.round(df.index + shift, 4)
                self.issues.append(Issue("major", "depth", f"Logging run {name} is off depth by {-shift:+.2f} {unit} "
                                         f"relative to {merged_name} (GR correlation {corr:.2f}); shifted before splicing",
                                         w.name, source=name))
                self.lineage.add("condition", "depth match run", w.name, gr, len(df), shift=round(shift, 3),
                                 corr=round(corr, 3), run=name)
            out_runs.append(df)
            merged_name = name
            merged = df
        spliced = splice(out_runs, priority="shallow_first")
        self.lineage.add("condition", "splice runs", w.name, n_affected=len(runs),
                         runs="|".join(r[0] for r in runs))
        self.issues.append(Issue("info", "format", f"{len(runs)} logging runs spliced", w.name))
        return spliced.drop(columns="RUN")

    @staticmethod
    def _name_from_filename(stem: str, names: list[str]) -> str | None:
        compact = re.sub(r"[\s_\-]", "", stem.upper())
        for n in sorted(names, key=len, reverse=True):
            key = re.sub(r"[\s_\-]", "", n.upper())
            root, num = re.match(r"^(.*?)(\d*)$", key).groups()
            pat = re.escape(root) + r"0*" + (num.lstrip("0") or "0") + r"(?!\d)" if num else re.escape(key)
            if re.search(pat, compact):
                return n
        return stem

    def _add(self, issues: list[Issue], proj: str, well: str | None = None) -> None:
        for i in issues:
            i.well = i.well or well
            i.source = f"{proj}/{i.source}" if i.source else proj
        self.issues.extend(issues)

    def _check_consistency(self, w: Well, meta: dict, source: str) -> None:
        if meta.get("kb") is not None and w.kb is not None and not pd.isna(w.kb) and abs(meta["kb"] - w.kb) > 0.1:
            self.issues.append(Issue("major", "datum", f"KB {meta['kb']:g} in deviation file vs {w.kb:g} in header",
                                     w.name, source=source))
        if None not in (meta.get("x"), w.x) and np.hypot(meta["x"] - w.x, meta["y"] - w.y) > 1.0:
            self.issues.append(Issue("major", "location", "Wellhead XY differs between deviation file and header "
                                     f"by {np.hypot(meta['x'] - w.x, meta['y'] - w.y):.1f}", w.name, source=source))
        du = meta.get("depth_unit")
        if du and w.depth_unit and du != w.depth_unit:
            self.issues.append(Issue("critical", "units", f"Deviation depth unit '{du}' vs '{w.depth_unit}'",
                                     w.name, source=source))
        w.depth_unit = w.depth_unit or du
        if meta.get("well_type"):
            w.sources.append(f"type={meta['well_type']}")

    # ------------------------------------------------------------------ standardise
    def standardise(self) -> None:
        for w in self.wells:
            if not w.depth_unit:
                self.issues.append(Issue("critical", "units", "Depth unit unknown; assumed target unit", w.name))
                w.depth_unit = self.depth_unit
            to_depth_unit(w, self.depth_unit, self.lineage)
            self.issues.extend(standardise_curves(w, self.lineage))
            if w.trajectory is not None and w.td is not None and not pd.isna(w.td):
                d = float(w.trajectory["MD"].max() - w.td)
                if abs(d) > 1.0:
                    self.issues.append(Issue("major", "depth", f"Deviation survey ends at {w.trajectory['MD'].max():.1f} "
                                             f"but header TD is {w.td:.1f} {w.depth_unit}", w.name))
        self.duplicates = find_duplicate_wells(self.wells)
        if not self.duplicates.empty:
            for _, r in self.duplicates[self.duplicates.verdict == "duplicate"].iterrows():
                self.issues.append(Issue("critical", "duplicate",
                                         f"Same physical well as {r.well_b} ({r.project_b}, UWI {r.uwi_b}): "
                                         f"identical location, TD and trajectory; curves identical "
                                         f"{r.curves_identical}", r.well_a))
        close = find_close_wells(self.wells)
        if not close.empty:
            dup_pairs = set()
            if not self.duplicates.empty:
                dup_pairs = {frozenset((r.well_a, r.well_b)) for _, r in self.duplicates.iterrows()
                             if r.verdict == "duplicate"}
            for _, r in close.iterrows():
                if frozenset((r.well_a, r.well_b)) not in dup_pairs and r.surface_dist > 5:
                    self.issues.append(Issue("info", "location", f"Wellhead {r.surface_dist:.0f} m from {r.well_b} "
                                                                 "(same pad?)", r.well_a))
        self.master = assign_master_uwis(self.wells, self.duplicates, self.uwi_prefix, self.lineage)

    # ------------------------------------------------------------------ condition
    def condition(self) -> pd.DataFrame:
        scores = []
        for w in self.wells:
            if w.logs is None:
                self.issues.append(Issue("major", "inventory", "No log data for well", w.name))
                continue
            clean, iss = condition_well(w, self.lineage)
            self.issues.extend(iss)
            self.raw_logs[w.name] = w.logs
            self.clean_logs[w.name] = clean
            scores.append(curve_scores(w.logs, clean, w.name))
        return pd.concat(scores, ignore_index=True) if scores else pd.DataFrame()

    # ------------------------------------------------------------------ checkshots & seismic
    def qc_time_depth(self) -> pd.DataFrame:
        rows = []
        for w in self.wells:
            if w.checkshot is None:
                continue
            df, iss = qc_checkshot(w.checkshot, w.name)
            self.issues.extend(iss)
            self.checkshot_qc[w.name] = df
            v = df["VINT"].dropna()
            rows.append({"well": w.name, "levels": len(df), "top": df.DEPTH.min(), "base": df.DEPTH.max(),
                         "vint_min": round(v.min()), "vint_max": round(v.max()),
                         "vint_smooth_range": f"{df.VINT_SMOOTH.min():.0f}-{df.VINT_SMOOTH.max():.0f}",
                         "flagged": int((df.FLAG != "").sum())})
        return pd.DataFrame(rows)

    def qc_seismic(self) -> pd.DataFrame:
        rows = []
        for proj, s in self.segy:
            rows.append({"project": proj, "file": Path(s.path).name, "inlines": f"{s.inline_range}",
                         "xlines": f"{s.xline_range}", "il/xl bytes": f"{s.inline_byte}/{s.xline_byte}",
                         "xy bytes": f"{s.x_byte}/{s.y_byte}",
                         "bin (m)": f"{s.grid['il_spacing']:.1f} x {s.grid['xl_spacing']:.1f}" if s.grid else "",
                         "time (ms)": f"{s.start_time_ms:g}-{s.start_time_ms + (s.n_samples - 1) * s.sample_interval_us / 1000:g}",
                         "traces": s.n_traces_complete, "missing": s.missing_traces})
        for i in range(len(self.segy)):
            for j in range(i + 1, len(self.segy)):
                a, b = self.segy[i][1], self.segy[j][1]
                if same_grid(a, b):
                    self.issues.append(Issue("info", "segy", f"{Path(a.path).name} and {Path(b.path).name} share one "
                                                             "inline/crossline grid (crops/versions of the same survey)"))
        # wells outside seismic coverage
        for w in self.wells:
            if w.x is None or not self.segy:
                continue
            x, y = (w.trajectory.X.iloc[-1], w.trajectory.Y.iloc[-1]) if w.trajectory is not None else (w.x, w.y)
            inside = []
            for _, s in self.segy:
                if not s.grid:
                    continue
                il, xl = xy_to_grid(s.grid, [x], [y])
                inside.append(s.inline_range[0] <= il[0] <= s.inline_range[1]
                              and s.xline_range[0] <= xl[0] <= s.xline_range[1])
            near = any(s.grid and s.grid["xmin"] - 10_000 < x < s.grid["xmax"] + 10_000
                       and s.grid["ymin"] - 10_000 < y < s.grid["ymax"] + 10_000 for _, s in self.segy)
            if inside and near and not any(inside):
                self.issues.append(Issue("major", "location", "Bottom-hole location outside all seismic surveys",
                                         w.name))
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------ prediction
    def predict(self, features=("GR", "RHOB", "NPHI", "TVDSS"), area_km: float = 20.0) -> None:
        """Predict DT where missing (per spatial area, one copy per physical well), then Vs everywhere."""
        reps: dict[str, Well] = {}
        for w in self.wells:
            if w.name not in self.clean_logs:
                continue
            k = w.master_uwi or w.name
            cur = reps.get(k)
            if cur is None or self.clean_logs[w.name].notna().sum().sum() > self.clean_logs[cur.name].notna().sum().sum():
                reps[k] = w
        areas = _cluster([w for w in reps.values() if w.x is not None], area_km * 1000)
        vals = []
        for area in areas:
            frames = {}
            for w in area:
                d = self.clean_logs[w.name].copy()
                d["TVDSS"] = w.md_to_tvdss(d.index.to_numpy())
                frames[w.name] = d
            try:
                preds, val, _, blind = predict_curve(frames, "DT", list(features))
                self.blind_predictions.update(blind)
            except ValueError as e:
                self.issues.append(Issue("info", "prediction", f"DT prediction skipped for area of "
                                                               f"{[w.name for w in area]}: {e}"))
                continue
            val["area_wells"] = len(area)
            small = val[val.n < 500]
            for _, r in small.iterrows():
                self.issues.append(Issue("info", "prediction", f"Blind test skipped: only {r.n} samples where DT "
                                                               "and all features overlap", r.blind_well, "DT"))
            val = val[val.n >= 500]
            for _, r in val[val.r2 < 0.3].iterrows():
                self.issues.append(Issue("major", "prediction",
                                         f"DT from neighbouring wells fails the blind test (R² {r.r2:.2f}, RMSE "
                                         f"{r.rmse:.1f} us/ft): this well's DT or input curves are inconsistent with "
                                         "the rest of the area; review scaling/calibration before use", r.blind_well, "DT"))
            vals.append(val)
            for w in area:
                d = frames[w.name]
                p = preds[w.name]
                dt_meas = d["DT"] if "DT" in d else pd.Series(np.nan, index=d.index)
                fill = dt_meas.isna() & p.notna()
                dt_final = dt_meas.where(~fill, p)
                self.predictions.setdefault(w.name, {})["DT"] = p
                n_fill = int(fill.sum())
                if n_fill:
                    self.lineage.add("predict", "fill missing DT with ML prediction", w.name, "DT", n_fill,
                                     features="|".join(features))
                # Synthetic shear via Greenberg-Castagna, Vsh from a GR index
                if "GR" in d and dt_final.notna().any():
                    gr = d["GR"]
                    lo, hi = gr.quantile(0.05), gr.quantile(0.95)
                    vsh = ((gr - lo) / (hi - lo)).clip(0, 1)
                    vp = dt_to_vp(dt_final)
                    vs = greenberg_castagna_vs(vp, {"sandstone": 1 - vsh + 1e-6, "shale": vsh + 1e-6})
                    dts = pd.Series(vp_to_dt(vs), index=d.index, name="DTS_SYN")
                    self.predictions[w.name]["DTS"] = dts
                    self.predictions[w.name]["DT_FINAL"] = dt_final
                    self.predictions[w.name]["VSH_GR"] = vsh
                    self.lineage.add("predict", "synthetic Vs (Greenberg-Castagna, Vsh from GR index)", w.name,
                                     "DTS", int(dts.notna().sum()))
                    if fill.any():
                        self.issues.append(Issue("info", "prediction", f"DT predicted over {n_fill} samples "
                                                                       "with no measured sonic", w.name, "DT"))
        self.validation = pd.concat(vals, ignore_index=True) if vals else pd.DataFrame()

    def rpt_qc(self, tol: float = 0.12) -> pd.DataFrame:
        tmpl = rpt_lines()
        rows = []
        for name, p in self.predictions.items():
            if "DTS" not in p or "RHOB" not in self.clean_logs[name]:
                continue
            vp, vs = dt_to_vp(p["DT_FINAL"]), dt_to_vp(p["DTS"])
            rho = self.clean_logs[name]["RHOB"]
            df = pd.DataFrame({"AI": vp * rho, "VpVs": vp / vs, "well": name}).dropna()
            if df.empty:
                continue
            df["RPT_DIST"] = rpt_distance(df.AI, df.VpVs, tmpl)
            df["RPT_FLAG"] = df.RPT_DIST > tol
            rows.append(df)
            n = int(df.RPT_FLAG.sum())
            if n / len(df) > 0.05:
                self.issues.append(Issue("major", "rock physics", f"{n} of {len(df)} samples ({n / len(df):.0%}) plot "
                                                                  "off the rock physics template; check DT/RHOB", name))
        return pd.concat(rows) if rows else pd.DataFrame()

    # ------------------------------------------------------------------ outputs
    def write(self, out: Path, scores: pd.DataFrame, cs_table: pd.DataFrame, seis_table: pd.DataFrame,
              rpt: pd.DataFrame, make_figures: bool = True) -> Path:
        out.mkdir(parents=True, exist_ok=True)
        (out / "las").mkdir(exist_ok=True)
        figs = out / "figures"
        figs.mkdir(exist_ok=True)
        inv = pd.DataFrame(self.inventory).drop(columns="path")
        issues = issues_frame(self.issues)
        lineage = self.lineage.to_frame()
        inv.to_csv(out / "inventory.csv", index=False)
        self.master.to_csv(out / "well_master.csv", index=False)
        self.duplicates.to_csv(out / "duplicates.csv", index=False)
        issues.to_csv(out / "issues.csv", index=False)
        scores.to_csv(out / "curve_scores.csv", index=False)
        lineage.to_csv(out / "lineage.csv", index=False)
        seis_table.to_csv(out / "segy.csv", index=False)
        if not self.validation.empty:
            self.validation.to_csv(out / "prediction_validation.csv", index=False)
        for name, df in self.checkshot_qc.items():
            df.to_csv(out / f"checkshot_qc_{canonical_name(name)}.csv", index=False)

        written = set()
        for w in self.wells:
            if w.name not in self.clean_logs or w.master_uwi in written:
                continue
            written.add(w.master_uwi)
            self._write_las(w, out / "las" / f"{w.master_uwi}_{canonical_name(w.name)}.las")

        fig_paths: list[str] = []
        if make_figures:
            fig_paths = self._figures(figs, scores, rpt)
        self.figures = fig_paths
        return write_report(out, inv, self.master, self.duplicates, issues, scores, seis_table, cs_table,
                            self.validation, lineage, fig_paths)

    def _write_las(self, w: Well, path: Path) -> None:
        las = lasio.LASFile()
        las.well["WELL"].value = canonical_name(w.name)
        las.well["UWI"].value = w.master_uwi or ""
        las.well["FLD"].value = w.project or ""
        las.well["COMP"].value = ""
        las.params.append(lasio.HeaderItem("ALIAS", value=";".join(sorted(w.aliases | {w.name})),
                                           descr="Legacy names"))
        las.params.append(lasio.HeaderItem("EKB", unit=w.depth_unit, value=w.kb if w.kb is not None else -999.25,
                                           descr="Kelly bushing elevation"))
        clean = self.clean_logs[w.name]
        las.append_curve("DEPT", clean.index.to_numpy(), unit=w.depth_unit)
        for c in clean.columns:
            las.append_curve(c, clean[c].to_numpy(), unit=w.units.get(c, ""), descr="conditioned")
        for c, s in self.predictions.get(w.name, {}).items():
            if c in ("DT", "DTS"):
                las.append_curve(f"{c}_PRED" if c == "DT" else "DTS_SYN", s.reindex(clean.index).to_numpy(),
                                 unit="us/ft", descr="ML prediction" if c == "DT" else "Greenberg-Castagna synthetic")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            las.write(str(path), version=2.0)

    def _figures(self, figs: Path, scores: pd.DataFrame, rpt: pd.DataFrame) -> list[str]:
        out = []
        if self.segy or not self.master.empty:
            surveys = [{"name": Path(s.path).stem, "grid": s.grid, "il": s.inline_range, "xl": s.xline_range}
                       for _, s in self.segy if s.grid]
            for area in _cluster([w for w in self.wells if w.x is not None], 20000):
                names = {w.name for w in area}
                m = self.master[self.master["aliases"].apply(lambda a, names=names: bool(set(a.split(";")) & names))]
                xs, ys = m.x.to_numpy(float), m.y.to_numpy(float)
                sv = [s for s in surveys if s["grid"]["xmin"] - 5000 < xs.mean() < s["grid"]["xmax"] + 5000
                      and s["grid"]["ymin"] - 5000 < ys.mean() < s["grid"]["ymax"] + 5000]
                tag = canonical_name(sorted(names)[0]).split("-")[0].lower()
                out.append(str(plots.basemap(m, sv, figs / f"basemap_{tag}.png",
                                             title=f"Wells and seismic coverage ({len(m)} physical wells)")))
        if not scores.empty and self.clean_logs:
            # showcase the well with the widest variety of conditioning fixes on key curves
            fixes = pd.Series([(i.well, i.category) for i in self.issues
                               if i.category in ("nulls", "spike", "flatline", "range")
                               and i.curve in ("GR", "DT", "RHOB", "NPHI") and i.well in self.clean_logs])
            if not fixes.empty:
                worst = fixes.drop_duplicates().map(lambda t: t[0]).value_counts().index[0]
            else:
                worst = next(iter(self.clean_logs))
            out.append(str(plots.log_panel(self.raw_logs[worst], self.clean_logs[worst], worst,
                                           self.depth_unit, figs / "log_conditioning.png",
                                           pred={"DT": self.predictions.get(worst, {}).get("DT")}
                                           if self.predictions.get(worst, {}).get("DT") is not None else None)))
        if not self.validation.empty:
            # show the median-scoring blind well: representative, not cherry-picked
            v = self.validation.sort_values("r2").reset_index(drop=True)
            pick = v.iloc[len(v) // 2]
            meas = self.clean_logs[pick.blind_well]["DT"]
            out.append(str(plots.blind_well(meas, self.blind_predictions[pick.blind_well].reindex(meas.index),
                                            pick.blind_well, pick.r2, self.depth_unit,
                                            figs / "dt_blind_well.png")))
        if not rpt.empty:
            out.append(str(plots.rpt_crossplot(rpt.sample(min(len(rpt), 40000), random_state=0), rpt_lines(),
                                               figs / "rpt_qc.png")))
        for name, df in list(self.checkshot_qc.items())[:2]:
            out.append(str(plots.checkshot_qc(df, name, self.depth_unit,
                                              figs / f"checkshot_{canonical_name(name)}.png")))
        return out


def _cluster(wells: list[Well], dist: float) -> list[list[Well]]:
    groups: list[list[Well]] = []
    for w in wells:
        hit = [g for g in groups if any(np.hypot(w.x - o.x, w.y - o.y) <= dist for o in g)]
        merged = [w] + [o for g in hit for o in g]
        groups = [g for g in groups if g not in hit] + [merged]
    return groups


def run_audit(projects: dict[str, str | Path], out_dir: str | Path, depth_unit: str = "m",
              predict: bool = True, figures: bool = True, uwi_prefix: str = "W") -> Audit:
    """Run the full audit and write all outputs to ``out_dir``. Returns the :class:`Audit` object."""
    a = Audit(projects, depth_unit, uwi_prefix)
    a.discover()
    a.read()
    a.standardise()
    scores = a.condition()
    cs = a.qc_time_depth()
    seis = a.qc_seismic()
    rpt = pd.DataFrame()
    if predict:
        a.predict()
        rpt = a.rpt_qc()
    a.report_path = a.write(Path(out_dir), scores, cs, seis, rpt, make_figures=figures)
    return a
