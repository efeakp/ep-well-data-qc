"""Synthetic legacy dataset with realistic defects, so the pipeline can be demonstrated and tested
without any proprietary data.

Two 'legacy projects' over one fictional field (NOVA):

* ``ALPHA_2009``: Petrel-style exports in **feet**, 4 wells
* ``BETA_2014``: exports in **metres**, 3 wells (one of which is ALPHA-03 reloaded under a new name/UWI),
  plus well tops, a checkshot and a small 3D SEG-Y cube

Injected defects (each one is something the audit is expected to find):

==  =========================================================================================
1   BETA 7 is the same physical well as ALPHA-03 (same location, TD, trajectory and curves)
2   ALPHA-02 DT contains sign-flipped nulls (+999.25); RHOB contains spikes
3   ALPHA-01 GR holds a constant filler value over 120 m
4   ALPHA-02 logs run 18 m below the header TD
5   Deviation files name wells 'ALPHAA-02' (typo) and 'CORE-04' (legacy name, location only)
6   ALPHA VSH labelled '%' but stored as a fraction; BETA NPHI labelled m3/m3 but stored in %
7   ALPHA-04 has no sonic (DT) at all; no well has shear (DTS)
8   BETA 6 was logged in two runs; the deep run is 2.4 m off depth
9   BETA checkshot has TIME/DEPTH columns reversed, 1 ms rounding and one mis-picked level
10  SEG-Y uses non-standard header bytes (IL @5, XL @21, XY @73/77) and is truncated
11  Deviation files state no CRS; ALPHA well header rows have extra columns
==  =========================================================================================
"""

from __future__ import annotations

import struct
from pathlib import Path

import lasio
import numpy as np
import pandas as pd

from .rockphysics import fluid_substitute, greenberg_castagna_vs, vp_to_dt

FT = 0.3048
STEP_M = 0.1524  # 0.5 ft

# Layer model: (top TVDSS m, base TVDSS m, facies)
LAYERS = [(0, 900, "shale"), (900, 980, "sand"), (980, 1350, "shale"), (1350, 1420, "sand"),
          (1420, 1700, "shale"), (1700, 1745, "sand"), (1745, 2050, "shale"), (2050, 2130, "gas sand"),
          (2130, 2400, "shale"), (2400, 2470, "sand"), (2470, 2800, "shale"), (2800, 2860, "sand"),
          (2860, 3300, "shale")]

WELLS = {
    # name: x, y, kb(m), td_md(m), max_incl(deg), azimuth, project
    "ALPHA-01": (512_340.0, 6_201_220.0, 25.0, 2950.0, 0.0, 0.0),
    "ALPHA-02": (513_910.0, 6_200_480.0, 26.0, 3050.0, 18.0, 120.0),
    "ALPHA-03": (511_620.0, 6_202_150.0, 25.0, 3000.0, 24.0, 300.0),
    "ALPHA-04": (514_480.0, 6_202_620.0, 24.0, 2900.0, 12.0, 45.0),
    "BETA 5": (512_980.0, 6_203_340.0, 18.0, 3100.0, 8.0, 200.0),
    "BETA 6": (513_350.0, 6_201_700.0, 18.0, 3080.0, 15.0, 80.0),
}


def _trajectory(x0, y0, kb, td, max_incl, azi, kop=600.0, build=2.0 / 30.0, step=30.0):
    md = np.arange(0.0, td + step, step)
    md[-1] = td
    incl = np.clip((md - kop) * build, 0, max_incl)
    a = np.radians(azi)
    dmd = np.diff(md, prepend=0)
    inc_mid = np.radians(np.r_[incl[0], (incl[1:] + incl[:-1]) / 2])
    dz = dmd * np.cos(inc_mid)
    dh = dmd * np.sin(inc_mid)
    tvd = np.cumsum(dz)
    dx, dy = np.cumsum(dh * np.sin(a)), np.cumsum(dh * np.cos(a))
    dls = np.r_[0, np.abs(np.diff(incl)) / np.diff(md) * 30]
    return pd.DataFrame({"MD": md, "X": x0 + dx, "Y": y0 + dy, "Z": kb - tvd, "TVD": tvd, "DX": dx, "DY": dy,
                         "AZIM": np.where(incl > 0, azi, 0.0), "INCL": incl, "DLS": dls})


def _earth(tvdss, x, y, rng):
    """Elastic and log properties at TVDSS (m) for one well; layers dip gently to the east."""
    shift = (x - 512_000) * 0.02 + (y - 6_200_000) * -0.01
    z = tvdss - shift
    n = z.size
    facies = np.full(n, "shale", dtype=object)
    for top, base, f in LAYERS:
        facies[(z >= top) & (z < base)] = f
    sand = facies != "shale"
    vsh = np.where(sand, 0.12, 0.82) + rng.normal(0, 0.05, n)
    vsh = np.clip(pd.Series(vsh).rolling(9, center=True, min_periods=1).mean().to_numpy(), 0.02, 0.98)
    zk = np.clip(z, 0, None) / 1000
    phi_sand = 0.36 * np.exp(-0.42 * zk)
    phi_sh = 0.42 * np.exp(-0.75 * zk)
    phi = np.where(sand, phi_sand, phi_sh) + rng.normal(0, 0.008, n)
    phi = np.clip(phi, 0.03, 0.4)
    # Vp: sands from Raymer-Hunt-Gardner, shales from a compaction trend
    vp_sand = (1 - phi) ** 2 * 5.6 + phi * 1.6
    vp_sh = 1.75 + 1.05 * zk - 0.9 * (phi_sh - phi)
    vp = np.where(sand, vp_sand, vp_sh)
    # geological variability the other logs cannot see: cementation/pressure effects that vary
    # smoothly with depth and from well to well, so predicted sonic is good but not perfect
    wobble = pd.Series(rng.normal(0, 1, n)).rolling(400, center=True, min_periods=1).mean().to_numpy()
    vp = vp * (1 + 0.05 * wobble / (np.std(wobble) or 1) + rng.normal(0, 0.02))
    vs = greenberg_castagna_vs(vp, {"sandstone": 1 - vsh, "shale": vsh})
    rho_ma = 2.65 * (1 - vsh) + 2.72 * vsh
    rho = (1 - phi) * rho_ma + phi * 1.05
    gas = facies == "gas sand"
    if gas.any():
        vp_g, vs_g, rho_g = fluid_substitute(vp[gas], vs[gas], rho[gas], phi[gas], 2.8, 1.05, 0.1, 0.2)
        vp[gas], vs[gas], rho[gas] = vp_g, vs_g, rho_g
    sw = np.where(gas, 0.25, 1.0)
    rw = 0.06
    rt = rw / (np.maximum(phi, 0.03) ** 2 * sw ** 2) * (1 + 2.5 * vsh)
    gr = 18 + 115 * vsh + rng.normal(0, 4, n)
    nphi = phi + 0.22 * vsh - np.where(gas, 0.12, 0) + rng.normal(0, 0.012, n)
    cali = 8.5 + 1.2 * vsh * (zk < 1.5) + np.abs(rng.normal(0, 0.15, n))
    dt = vp_to_dt(vp) + rng.normal(0, 2.0, n)
    dts = vp_to_dt(vs) + rng.normal(0, 2.5, n)
    rhob = rho + rng.normal(0, 0.012, n)
    return pd.DataFrame({"GR": gr, "DT": dt, "DTS": dts, "RHOB": rhob, "NPHI": nphi,
                         "ILD": rt * np.exp(rng.normal(0, 0.08, n)), "CALI": cali, "VSH": vsh})


def _write_las(path: Path, well: str, uwi: str, depth, df: pd.DataFrame, units: dict, depth_unit: str):
    las = lasio.LASFile()
    las.well["WELL"].value = well
    las.well["UWI"].value = uwi
    las.well["NULL"].value = -999.25
    las.append_curve("DEPT", np.round(depth, 4), unit=depth_unit)
    for c in df.columns:
        las.append_curve(c, df[c].fillna(-999.25).to_numpy(), unit=units.get(c, ""))
    las.write(str(path), version=2.0, fmt="%.4f")


def _write_dev(path: Path, name: str, traj: pd.DataFrame, x, y, kb, unit: str, crs="<Undefined>"):
    lines = ["# WELL TRACE FROM PETREL", f"# WELL NAME:              {name}",
             f"# WELL HEAD X-COORDINATE: {x:.8f} (m)", f"# WELL HEAD Y-COORDINATE: {y:.8f} (m)",
             f"# WELL KB:                {kb:.8f} ({unit})", "# WELL TYPE:              OIL",
             "# MD AND TVD ARE REFERENCED (=0) AT KB AND INCREASE DOWNWARDS", "# ANGLES ARE GIVEN IN DEGREES",
             f"# XYZ TRACE IS GIVEN IN COORDINATE SYSTEM {crs}", "# AZIMUTH REFERENCE UNDEFINED",
             "# DX DY ARE GIVEN IN GRID NORTH IN m-UNITS", f"# DEPTH (Z, TVD) GIVEN IN {unit}-UNITS",
             "#" + "=" * 100,
             "      MD              X              Y             Z           TVD           DX           DY"
             "          AZIM          INCL          DLS", "#" + "=" * 100]
    for _, r in traj.iterrows():
        lines.append(" ".join(f"{v:14.6f}" for v in r.to_numpy()))
    path.write_text("\n".join(lines) + "\n")


def _write_segy(path: Path, rng, il0=1001, xl0=2001, n_il=56, n_xl=46, ns=376, dt_us=4000, t0=1500,
                origin=(510_000.0, 6_199_000.0), bin_m=(100.0, 100.0), truncate_traces=37):
    text = "".join(f"{line:<80}" for line in [
        "C 1 SYNTHETIC 3D SEISMIC  NOVA_PSTM  (generated by wellqc.synthetic)",
        f"C 2 First inline: {il0}   Last inline: {il0 + n_il - 1}",
        f"C 3 First xline:  {xl0}   Last xline:  {xl0 + n_xl - 1}",
        "C 4 Inline bytes 5-8, Xline bytes 21-24, X/Y bytes 73/77, scalar -10"] + [f"C{i:2d}" for i in range(5, 41)])
    binh = bytearray(400)
    struct.pack_into(">h", binh, 16, dt_us)
    struct.pack_into(">h", binh, 20, ns)
    struct.pack_into(">h", binh, 24, 5)
    t = np.arange(ns) * dt_us / 1e6
    w = (1 - 2 * (np.pi * 30 * (t - 0.05)) ** 2) * np.exp(-((np.pi * 30 * (t - 0.05)) ** 2))
    refl = rng.normal(0, 1, ns) * (rng.random(ns) < 0.08)
    with open(path, "wb") as fh:
        fh.write(text.encode("cp500"))
        fh.write(binh)
        n_total = n_il * n_xl
        for k in range(n_total - truncate_traces):
            i, j = divmod(k, n_xl)
            th = bytearray(240)
            struct.pack_into(">i", th, 4, il0 + i)
            struct.pack_into(">i", th, 20, xl0 + j)
            x = origin[0] + i * bin_m[0]
            y = origin[1] + j * bin_m[1]
            struct.pack_into(">h", th, 70, -10)
            struct.pack_into(">i", th, 72, int(round(x * 10)))
            struct.pack_into(">i", th, 76, int(round(y * 10)))
            struct.pack_into(">h", th, 108, t0)
            struct.pack_into(">h", th, 114, ns)
            struct.pack_into(">h", th, 116, dt_us)
            r = np.roll(refl, int(3 * np.sin(i / 7) + 2 * np.cos(j / 5)))
            tr = np.convolve(r, w[:40], mode="same") * 1000 + rng.normal(0, 20, ns)
            fh.write(th)
            fh.write(tr.astype(">f4").tobytes())
        fh.write(bytes(500))  # partial trace at the end: the file was cut off mid-transfer


def make_demo_dataset(out: str | Path, seed: int = 7) -> dict[str, Path]:
    """Write the synthetic legacy projects to ``out`` and return {project name: folder}."""
    out = Path(out)
    rng = np.random.default_rng(seed)
    pa, pb = out / "ALPHA_2009", out / "BETA_2014"
    for p in (pa / "logs", pa / "deviation", pb / "las", pb / "dev", pb / "seismic"):
        p.mkdir(parents=True, exist_ok=True)

    logs_m, trajs = {}, {}
    for name, (x, y, kb, td, inc, azi) in WELLS.items():
        traj = _trajectory(x, y, kb, td, inc, azi)
        md = np.round(np.arange(round(300 / STEP_M) * STEP_M, td + STEP_M / 2, STEP_M), 4)
        tvd = np.interp(md, traj.MD, traj.TVD)
        xs, ys = np.interp(md, traj.MD, traj.X), np.interp(md, traj.MD, traj.Y)
        df = _earth(tvd - kb, xs, ys, rng)
        df.index = md
        logs_m[name], trajs[name] = df, traj

    # ---------------- ALPHA (feet)
    alpha = [n for n in WELLS if n.startswith("ALPHA")]
    head = ["# Petrel well head", "# Unit in X and Y direction: m", "# Unit in depth: ft", "VERSION 1",
            "BEGIN HEADER", "Name", "UWI", "Well symbol", "Surface X", "Surface Y", "KB", "TD (MD)", "Operator",
            "END HEADER"]
    for k, n in enumerate(alpha, 1):
        x, y, kb, td, *_ = WELLS[n]
        td_hdr = td - 18 if n == "ALPHA-02" else td  # defect 4: header TD shallower than logs
        head.append(f'"{n}" "NV{k:04d}-A" 3 {x:.2f} {y:.2f} {kb / FT:.2f} {td_hdr / FT:.1f} "" -999 5 0')
    (pa / "well_header.txt").write_text("\n".join(head) + "\n")  # defect 11: extra trailing values

    for n in alpha:
        x, y, kb, td, *_ = WELLS[n]
        df = logs_m[n].drop(columns=["DTS"]).copy()
        if n == "ALPHA-01":  # defect 3: filler GR
            sl = (df.index > 1000) & (df.index < 1120)
            df.loc[sl, "GR"] = 77.94
        if n == "ALPHA-02":  # defect 2: +999.25 and spikes
            idx = rng.choice(len(df), 25, replace=False)
            df.iloc[idx, df.columns.get_loc("DT")] = 999.25
            idx = rng.choice(len(df), 40, replace=False)
            df.iloc[idx, df.columns.get_loc("RHOB")] += rng.choice([-0.9, 0.7], 40)
        if n == "ALPHA-04":  # defect 7: no sonic
            df = df.drop(columns=["DT"])
        df = df.rename(columns={"VSH": "VSH_CED"})
        units = {"GR": "gAPI", "DT": "us/ft", "RHOB": "g/cm3", "NPHI": "v/v", "ILD": "ohm.m", "CALI": "in",
                 "VSH_CED": "%"}  # defect 6: % label on a fraction
        _write_las(pa / "logs" / f"{n}.las", n, f"NV{alpha.index(n) + 1:04d}-A", df.index / FT, df, units, "ft")
        dev_name = {"ALPHA-02": "ALPHAA-02", "ALPHA-04": "CORE-04"}.get(n, n)  # defect 5
        t = trajs[n].copy()
        for c in ("MD", "Z", "TVD"):
            t[c] = t[c] / FT
        _write_dev(pa / "deviation" / f"{n}.dev", dev_name, t, x, y, kb / FT, "ft")

    # ---------------- BETA (metres): BETA 5, BETA 6 and BETA 7 (= ALPHA-03)
    beta = {"BETA 5": "BETA 5", "BETA 6": "BETA 6", "BETA 7": "ALPHA-03"}
    rows = ["Name\tUWI\tSurface\tSurface\tKB\tTD"]
    for k, (n, src) in enumerate(beta.items(), 5):
        x, y, kb, td, *_ = WELLS[src]
        rows.append(f"{n}\tNVB{k:02d}\t{x:.2f}\t{y:.2f}\t{kb:g}\t{td:g}")
    (pb / "well_header.txt").write_text("\n".join(rows) + "\n")
    for k, (n, src) in enumerate(beta.items(), 5):
        x, y, kb, td, *_ = WELLS[src]
        df = logs_m[src].drop(columns=["VSH", "DTS"]).copy()
        df["NPHI"] = df["NPHI"] * 100  # defect 6: percent stored under a fraction label
        units = {"GR": "gAPI", "DT": "us/ft", "RHOB": "g/cm3", "NPHI": "m3/m3", "ILD": "ohm.m", "CALI": "in"}
        if n == "BETA 6":  # defect 8: two runs, deep run off depth by +2.4 m
            split = 1900.0
            r1 = df[df.index <= split + 60]
            r2 = df[df.index >= split - 60].copy()
            r2.index = np.round(r2.index + 2.4, 4)
            _write_las(pb / "las" / "BETA 6_run1.las", n, f"NVB{k:02d}", r1.index, r1, units, "m")
            _write_las(pb / "las" / "BETA 6_run2.las", n, f"NVB{k:02d}", r2.index, r2, units, "m")
        else:
            _write_las(pb / "las" / f"{n}.las", n, f"NVB{k:02d}", df.index, df, units, "m")
        _write_dev(pb / "dev" / f"{n}.dev", n, trajs[src], x, y, kb, "m", crs="<Null>")

    # tops (Petrel format) picked on the layer model
    tops = ["# Petrel well tops", "# Unit in X and Y direction: m", "# Unit in depth: m", "VERSION 2",
            "BEGIN HEADER", "X", "Y", "Z", "TWT picked", "MD", "Type", "Surface", "Well", "Interpreter",
            "END HEADER"]
    for n, src in beta.items():
        t = trajs[src]
        kb = WELLS[src][2]
        for surf, z_top in (("Reservoir A Top", 2050), ("Reservoir A Base", 2130)):
            md_grid = np.linspace(0, t.MD.max(), 5000)
            xs, ys = np.interp(md_grid, t.MD, t.X), np.interp(md_grid, t.MD, t.Y)
            tvdss = np.interp(md_grid, t.MD, t.TVD) - kb
            zz = tvdss - ((xs - 512_000) * 0.02 + (ys - 6_200_000) * -0.01)
            i = int(np.argmin(np.abs(zz - z_top)))
            tops.append(f'{xs[i]:.2f} {ys[i]:.2f} {-tvdss[i]:.2f} -999 {md_grid[i]:.2f} Horizon "{surf}" "{n}" "demo"')
    (pb / "reservoir_tops.dat").write_text("\n".join(tops) + "\n")

    # checkshot for BETA 5: TIME first, rounded to 1 ms, one mis-picked level (defect 9)
    src = "BETA 5"
    d = logs_m[src]
    md = d.index.to_numpy()
    kb = WELLS[src][2]
    tvd = np.interp(md, trajs[src].MD, trajs[src].TVD)
    owt = np.cumsum(np.r_[0, np.diff(tvd)] / (304.8 / d["DT"].to_numpy()) / 1000)  # s, from first sample
    v_over = 1.75  # km/s above the logged interval
    owt += (tvd[0] - kb) / v_over / 1000
    levels = np.arange(500.0, 3050.0, 50.0)
    twt = np.round(2000 * np.interp(levels, md, owt))
    twt[24] += 9  # mis-picked level
    cs = ["2-way time\tMeasured Depth"] + [f"{t:.0f}\t{z:.0f}" for t, z in zip(twt, levels, strict=True)]
    (pb / "BETA 5 CS.txt").write_text("\n".join(cs) + "\n")

    _write_segy(pb / "seismic" / "NOVA_PSTM.segy", rng)  # defect 10
    return {"ALPHA_2009": pa, "BETA_2014": pb}
