"""LAS reading with clean-up of redundant index columns and header sanity checks."""

from __future__ import annotations

import warnings
from pathlib import Path

import lasio
import numpy as np
import pandas as pd

from ..model import Issue

FT_PER_M = 3.280839895


def read_las(path: str | Path) -> tuple[dict, pd.DataFrame, dict[str, str], list[Issue]]:
    """Read a LAS 2.0 file.

    Returns ``(header, logs, units, issues)`` where ``logs`` is indexed by depth and has
    redundant depth/index copies removed (Petrel exports often repeat DEPT two or three times,
    sometimes in a different unit).
    """
    path = Path(path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        las = lasio.read(str(path), ignore_header_errors=True)
    issues: list[Issue] = []
    w = las.well

    def hv(key):
        try:
            v = w[key].value
        except KeyError:
            return None
        return v if v not in ("", None) else None

    header = {
        "well": hv("WELL"), "uwi": hv("UWI"), "field": hv("FLD"), "company": hv("COMP"),
        "date": hv("DATE"), "strt": hv("STRT"), "stop": hv("STOP"), "step": hv("STEP"),
        "null": hv("NULL"), "depth_unit": (las.curves[0].unit or las.well["STRT"].unit or "").lower() or None,
        "source": path.name,
    }
    name = header["well"] or path.stem
    for k in ("field", "company"):
        if header[k] is None:
            issues.append(Issue("minor", "metadata", f"LAS header {k.upper()} is empty", well=name,
                                source=path.name))

    df = las.df()
    df.index = df.index.astype(float)
    df.index.name = "DEPT"
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
    # Clean the lasio duplicate suffix (e.g. 'DT:1') when it is now unique
    rename = {}
    for c in df.columns:
        base = c.split(":")[0]
        if ":" in c and sum(x.split(":")[0] == base for x in df.columns) == 1:
            rename[c] = base
    df = df.rename(columns=rename)
    units = {rename.get(k, k): v for k, v in units.items()}
    return header, df, units, issues
