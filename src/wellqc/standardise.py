"""Standardisation: curve mnemonics, units, well names, duplicate wells and master UWIs."""

from __future__ import annotations

import re
from difflib import SequenceMatcher
from itertools import combinations

import numpy as np
import pandas as pd

from .model import Issue, Lineage, Well, normalise_depth_unit

FT_PER_M = 3.280839895

# Canonical mnemonic -> accepted aliases (first match wins; later matches kept as <CANON>_ALT<n>)
CURVE_ALIASES: dict[str, list[str]] = {
    "GR": ["GR", "GRC", "GR_EDTC", "SGR_TOTAL"],
    "DT": ["DT", "DTC", "DTCO", "AC", "DTLF", "DT24"],
    "DTS": ["DTS", "DTSM", "DTSH"],
    "RHOB": ["RHOB", "RHOZ", "DEN", "ZDEN"],
    "NPHI": ["NPHI", "TNPH", "NPOR", "CNC"],
    "DRHO": ["DRHO", "HDRA"],
    "CALI": ["CALI", "CAL", "HCAL", "C1"],
    "RDEEP": ["ILD", "RILD", "LLD", "RT", "RDEP", "RD", "AT90", "HDRS"],
    "RMED": ["ILM", "LLS", "RILM", "RMED", "RM", "HMRS"],
    "RSHAL": ["MSFL", "SFLU", "SFL", "LLS3", "RXO", "SN", "RS", "RXOZ"],
    "SP": ["SP"],
    "PHIT": ["PHIT", "PHIE"],
    "PE": ["PE", "PEF", "PEFZ"],
    "VSH": ["VSH", "VSH_CED", "VCL", "VSHALE"],
    "URAN": ["URAN"],
    "THOR": ["THOR"],
    "POTA": ["POTA"],
    "SGR": ["SGR"],
}

# Physically plausible ranges in canonical units (values outside are set to null and counted)
PHYSICAL_RANGES: dict[str, tuple[float, float]] = {
    "GR": (0, 400), "PE": (0, 15), "DT": (40, 200), "DTS": (60, 650), "RHOB": (1.2, 3.2), "NPHI": (-0.15, 1.0),
    "DRHO": (-0.5, 0.5), "CALI": (4, 30), "RDEEP": (0.05, 5000), "RMED": (0.05, 5000),
    "RSHAL": (0.05, 5000), "SP": (-300, 300), "PHIT": (0, 0.6), "VSH": (0, 1),
    "URAN": (0, 100), "THOR": (0, 100), "POTA": (0, 10), "SGR": (0, 400),
}

CANONICAL_UNITS = {"GR": "gAPI", "PE": "b/e", "DT": "us/ft", "DTS": "us/ft", "RHOB": "g/cm3", "NPHI": "v/v",
                   "DRHO": "g/cm3", "CALI": "in", "RDEEP": "ohm.m", "RMED": "ohm.m", "RSHAL": "ohm.m",
                   "SP": "mV", "PHIT": "v/v", "VSH": "v/v", "URAN": "ppm", "THOR": "ppm", "POTA": "%",
                   "SGR": "gAPI"}

_FRACTION_CURVES = ("NPHI", "PHIT", "VSH")


# --------------------------------------------------------------------------- curves and units
def standardise_curves(well: Well, lineage: Lineage) -> list[Issue]:
    """Rename curves to canonical mnemonics and fix fraction/percent scaling."""
    issues: list[Issue] = []
    if well.logs is None:
        return issues
    lookup = {}
    for canon, aliases in CURVE_ALIASES.items():
        for a in aliases:
            lookup.setdefault(a.upper(), canon)
    rename, used = {}, {}
    for c in well.logs.columns:
        base = c.split(":")[0].upper()
        canon = lookup.get(base)
        if canon is None:
            continue
        if canon in used:
            used[canon] += 1
            rename[c] = f"{canon}_ALT{used[canon]}"
        else:
            used[canon] = 0
            rename[c] = canon
    well.logs = well.logs.rename(columns=rename)
    well.units = {rename.get(k, k): v for k, v in well.units.items()}
    if rename:
        lineage.add("standardise", "rename curves", well.name,
                    n_affected=len(rename), mapping="; ".join(f"{k}->{v}" for k, v in rename.items() if k != v))

    for c in list(well.logs.columns):
        canon = c.split("_ALT")[0]
        if canon not in _FRACTION_CURVES:
            continue
        s = well.logs[c].dropna()
        if s.empty:
            continue
        unit = (well.units.get(c) or "").lower()
        p99 = s.quantile(0.99)
        says_percent = unit in ("%", "pu", "percent")
        is_percent = p99 > 1.5
        if is_percent:
            well.logs[c] = well.logs[c] / 100.0
            lineage.add("standardise", "percent -> fraction", well.name, c, int(s.size), p99=round(p99, 2))
            if not says_percent:
                issues.append(Issue("major", "units",
                                    f"{c} labelled '{well.units.get(c)}' but values are in percent "
                                    f"(p99 = {p99:.1f}); divided by 100", well.name, c))
        elif says_percent:
            issues.append(Issue("minor", "units",
                                f"{c} labelled '%' but values are fractions (p99 = {p99:.2f}); label corrected",
                                well.name, c))
        well.units[c] = "v/v"
    return issues


def to_depth_unit(well: Well, target: str, lineage: Lineage) -> None:
    """Convert MD index, KB, TD, trajectory depths and tops to ``target`` ('m' or 'ft')."""
    src = normalise_depth_unit(well.depth_unit or target)
    target = normalise_depth_unit(target)
    if src is None or target is None:
        raise ValueError(f"{well.name}: cannot convert depth unit '{well.depth_unit}' to '{target}'")
    if src == target:
        well.depth_unit = target
        return
    f = 1 / FT_PER_M if (src, target) == ("ft", "m") else FT_PER_M
    if well.logs is not None:
        well.logs.index = np.round(well.logs.index.to_numpy() * f, 4)
    for attr in ("kb", "td"):
        v = getattr(well, attr)
        if v is not None and not pd.isna(v):
            setattr(well, attr, v * f)
    if well.trajectory is not None:
        for c in ("MD", "Z", "TVD"):
            if c in well.trajectory:
                well.trajectory[c] = well.trajectory[c] * f
    if well.tops is not None:
        for c in ("MD", "Z"):
            if c in well.tops:
                well.tops[c] = well.tops[c] * f
    well.depth_unit = target
    lineage.add("standardise", f"depth unit {src} -> {target}", well.name, factor=round(f, 6))


# --------------------------------------------------------------------------- well names
def name_key(name: str) -> tuple[str, int | None, str]:
    """Split a well name into (alpha root, number, suffix): 'DELTA-02ST1' -> ('DELTA', 2, 'ST1')."""
    s = re.sub(r"[\s_\-]+", " ", str(name).upper()).strip()
    m = re.match(r"^([A-Z ]*?)\s*0*(\d+)\s*([A-Z0-9 ]*)$", s)
    if not m:
        return s, None, ""
    return m.group(1).strip(), int(m.group(2)), m.group(3).replace(" ", "")


def canonical_name(name: str, width: int = 2) -> str:
    root, num, suf = name_key(name)
    return f"{root}-{num:0{width}d}{suf}" if num is not None else root


def match_name(name: str, candidates: list[str], min_ratio: float = 0.8) -> tuple[str | None, str]:
    """Match a well name to a list of known names.

    Returns (match, method) where method is 'exact', 'fuzzy' or 'none'. Fuzzy matching requires
    the same well number and a similar root, so 'ALPHAA 2' matches 'ALPHA 2' but not 'ALPHA 12'.
    """
    root, num, suf = name_key(name)
    for c in candidates:
        if name_key(c) == (root, num, suf):
            return c, "exact"
    best, best_r = None, 0.0
    for c in candidates:
        r2, n2, s2 = name_key(c)
        if n2 != num or s2 != suf:
            continue
        r = SequenceMatcher(None, root, r2).ratio()
        if r > best_r:
            best, best_r = c, r
    return (best, "fuzzy") if best_r >= min_ratio else (None, "none")


def match_by_location(x: float, y: float, heads: pd.DataFrame, tol: float = 2.0) -> str | None:
    d = np.hypot(heads["x"].to_numpy(float) - x, heads["y"].to_numpy(float) - y)
    i = int(np.argmin(d))
    return heads["name"].iloc[i] if d[i] <= tol else None


# --------------------------------------------------------------------------- duplicates
def find_duplicate_wells(wells: list[Well], xy_tol: float = 5.0, kb_tol: float = 0.5,
                         td_tol: float = 1.0) -> pd.DataFrame:
    """Find pairs of wells that are the same physical well loaded under different names/UWIs.

    Evidence: surface location, KB, TD, trajectory end point and identical curve samples.
    """
    rows = []
    for a, b in combinations(wells, 2):
        if None in (a.x, a.y, b.x, b.y):
            continue
        d = float(np.hypot(a.x - b.x, a.y - b.y))
        if d > xy_tol:
            continue
        ev = {"surface_dist": round(d, 2)}
        ev["kb_match"] = a.kb is not None and b.kb is not None and abs(a.kb - b.kb) <= kb_tol
        ev["td_match"] = a.td is not None and b.td is not None and abs(a.td - b.td) <= td_tol
        if a.trajectory is not None and b.trajectory is not None:
            ea, eb = a.trajectory.iloc[-1], b.trajectory.iloc[-1]
            ev["bottomhole_dist"] = round(float(np.hypot(ea.X - eb.X, ea.Y - eb.Y)), 2)
        identical, compared = [], []
        if a.logs is not None and b.logs is not None:
            j = a.logs.join(b.logs, how="inner", lsuffix="_a", rsuffix="_b")
            for c in set(a.logs.columns) & set(b.logs.columns):
                va, vb = j[f"{c}_a"], j[f"{c}_b"]
                ok = va.notna() & vb.notna()
                if ok.sum() < 50:
                    continue
                compared.append(c)
                if np.allclose(va[ok], vb[ok], rtol=1e-6, atol=1e-6):
                    identical.append(c)
        ev["curves_identical"] = f"{len(identical)}/{len(compared)}"
        same_hole = ev.get("bottomhole_dist", 0) < 1.0 and ev["td_match"]
        verdict = ("duplicate" if same_hole and (not compared or len(identical) >= max(1, len(compared) // 2))
                   else "same pad / sidetrack?" if d < xy_tol else "")
        rows.append({"well_a": a.name, "uwi_a": a.uwi, "project_a": a.project,
                     "well_b": b.name, "uwi_b": b.uwi, "project_b": b.project, **ev,
                     "identical_curves": ",".join(sorted(identical)), "verdict": verdict})
    return pd.DataFrame(rows)


def find_close_wells(wells: list[Well], max_dist: float = 100.0) -> pd.DataFrame:
    rows = []
    for a, b in combinations(wells, 2):
        if None in (a.x, a.y, b.x, b.y):
            continue
        d = float(np.hypot(a.x - b.x, a.y - b.y))
        if d <= max_dist:
            rows.append({"well_a": a.name, "well_b": b.name, "surface_dist": round(d, 1)})
    return pd.DataFrame(rows)


def assign_master_uwis(wells: list[Well], duplicates: pd.DataFrame, prefix: str = "W",
                       lineage: Lineage | None = None) -> pd.DataFrame:
    """Give every physical well one master UWI and keep all project names/UWIs as aliases."""
    parent = {w.name: w.name for w in wells}

    def find(n):
        while parent[n] != n:
            parent[n] = parent[parent[n]]
            n = parent[n]
        return n

    if not duplicates.empty:
        for _, r in duplicates[duplicates["verdict"] == "duplicate"].iterrows():
            parent[find(r["well_b"])] = find(r["well_a"])
    groups: dict[str, list[Well]] = {}
    for w in wells:
        groups.setdefault(find(w.name), []).append(w)
    rows = []
    for i, (_, members) in enumerate(sorted(groups.items(), key=lambda kv: (kv[1][0].project or "", kv[0])), 1):
        uwi = f"{prefix}{i:04d}"
        names = sorted({m.name for m in members} | {a for m in members for a in m.aliases})
        uwis = sorted({m.uwi for m in members if m.uwi})
        ref = members[0]
        for m in members:
            m.master_uwi = uwi
        rnd = lambda v, n=2: None if v is None or pd.isna(v) else round(float(v), n)  # noqa: E731
        rows.append({"master_uwi": uwi, "standard_name": canonical_name(ref.name), "x": rnd(ref.x), "y": rnd(ref.y),
                     "kb": rnd(ref.kb), "td": rnd(ref.td), "depth_unit": ref.depth_unit,
                     "projects": ";".join(sorted({m.project or "" for m in members})),
                     "aliases": ";".join(names), "legacy_uwis": ";".join(uwis), "n_sources": len(members)})
        if lineage is not None and len(members) > 1:
            lineage.add("standardise", "merge duplicate wells", ref.name, n_affected=len(members),
                        master_uwi=uwi, members="|".join(m.name for m in members))
    return pd.DataFrame(rows)
