"""LAS reading with clean-up of redundant index columns and header sanity checks."""

from __future__ import annotations

import warnings
from pathlib import Path

import lasio
import numpy as np
import pandas as pd

from ..model import Issue, normalise_depth_unit

FT_PER_M = 3.280839895

# Header mnemonics that carry well metadata, in order of preference
_X_KEYS = ("X", "XCOORD", "X_COORD", "XWELL", "EAST", "EASTING")
_Y_KEYS = ("Y", "YCOORD", "Y_COORD", "YWELL", "NORTH", "NORTHING")
_KB_KEYS = ("EKB", "KB", "EDF", "APD")
_TD_KEYS = ("TDD", "TDL", "TD")


def read_las(path: str | Path) -> tuple[dict, pd.DataFrame, dict[str, str], list[Issue]]:
    """Read a LAS 2.0 file.

    Returns ``(header, logs, units, issues)`` where ``logs`` is indexed by depth (always increasing)
    and has redundant depth/index copies removed (Petrel and DLIS-converted exports often repeat
    DEPT, sometimes in a different unit). ``header`` also carries any location, KB and TD found in
    the ~Well or ~Parameter sections, so a LAS file can stand in for a missing well header.
    """
    path = Path(path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        las = lasio.read(str(path), ignore_header_errors=True)
    issues: list[Issue] = []

    def item(keys):
        for section in (las.well, las.params):
            for k in keys:
                try:
                    it = section[k]
                except KeyError:
                    continue
                if it.value not in ("", None):
                    return it
        return None

    def hv(key):
        it = item((key,))
        return it.value if it is not None else None

    raw_unit = las.curves[0].unit or las.well["STRT"].unit or ""
    depth_unit = normalise_depth_unit(raw_unit)
    header = {
        "well": hv("WELL"), "uwi": hv("UWI"), "field": hv("FLD"), "company": hv("COMP"),
        "date": hv("DATE"), "strt": hv("STRT"), "stop": hv("STOP"), "step": hv("STEP"),
        "null": hv("NULL"), "depth_unit": depth_unit, "source": path.name,
    }
    name = header["well"] or path.stem
    if depth_unit is None:
        issues.append(Issue("critical", "units", f"Depth unit '{raw_unit}' not recognised", well=name,
                            source=path.name))
    if header["well"] is None:
        issues.append(Issue("major", "naming", f"LAS WELL field is empty; well named '{path.stem}' from the file name",
                            well=name, source=path.name))
    for k in ("field", "company"):
        if header[k] is None or str(header[k]).strip().upper() in ("UNKNOWN", "N/A", "NA"):
            issues.append(Issue("minor", "metadata", f"LAS header {k.upper()} is empty", well=name,
                                source=path.name))

    # Location / KB / TD from the header, converted to the log depth unit where they differ
    def num(keys, convert=True):
        it = item(keys)
        if it is None:
            return None, None
        try:
            v = float(it.value)
        except (TypeError, ValueError):
            return None, None
        u = normalise_depth_unit(it.unit)
        if convert and u and depth_unit and u != depth_unit:
            v = v * (FT_PER_M if (u, depth_unit) == ("m", "ft") else 1 / FT_PER_M)
        return v, it.mnemonic

    header["x"], _ = num(_X_KEYS, convert=False)
    header["y"], _ = num(_Y_KEYS, convert=False)
    header["lat"], _ = num(("LATI", "LAT"), convert=False)
    header["lon"], _ = num(("LONG", "LON"), convert=False)
    header["kb"], header["kb_key"] = num(_KB_KEYS)
    header["td"], header["td_key"] = num(_TD_KEYS)

    df = las.df()
    df.index = df.index.astype(float)
    df.index.name = "DEPT"
    if len(df) > 1 and df.index[0] > df.index[-1]:
        df = df.sort_index()
        issues.append(Issue("info", "format", "Recorded bottom-up (negative STEP); re-sorted to increasing depth",
                            well=name, source=path.name))
    units = {c.mnemonic: (c.unit or "") for c in las.curves[1:]}
    dropped = []
    idx = df.index.to_numpy()
    for c in list(df.columns):
        v = df[c].to_numpy(dtype=float)
        ok = ~np.isnan(v)
        if ok.sum() < 2:
            continue
        for factor, label in ((1.0, "copy of depth index"), (FT_PER_M, "depth index in ft"),
                              (1 / FT_PER_M, "depth index in m")):
            if np.allclose(v[ok], idx[ok] * factor, rtol=1e-6, atol=1e-3):
                dropped.append(f"{c} ({label})")
                df = df.drop(columns=c)
                units.pop(c, None)
                break
    if dropped:
        issues.append(Issue("info", "format", "Dropped redundant index curves: " + "; ".join(dropped),
                            well=name, source=path.name))
    empty = [c for c in df.columns if df[c].notna().sum() == 0]
    if empty:
        df = df.drop(columns=empty)
        issues.append(Issue("minor", "inventory", "Curves with no data dropped: " + ", ".join(empty),
                            well=name, source=path.name))
    # Clean the lasio duplicate suffix (e.g. 'DT:1') when it is now unique
    rename = {}
    for c in df.columns:
        base = c.split(":")[0]
        if ":" in c and sum(x.split(":")[0] == base for x in df.columns) == 1:
            rename[c] = base
    df = df.rename(columns=rename)
    units = {rename.get(k, k): v for k, v in units.items() if rename.get(k, k) in df.columns}
    return header, df, units, issues
