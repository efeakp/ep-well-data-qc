"""Log conditioning: nulls, physical limits, despiking, flat-line detection, depth matching, splicing.

Every function is non-destructive (returns a new Series/DataFrame) and reports what it changed so
the caller can record it in the :class:`~wellqc.model.Lineage`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .model import Issue, Lineage, Well
from .standardise import PHYSICAL_RANGES

NULL_SENTINELS = (-999.25, 999.25, -999.0, -9999.0, -99999.0, 1e30, -1e30)


def fix_nulls(s: pd.Series) -> tuple[pd.Series, int]:
    """Replace null sentinels that slipped through (including the sign-flipped +999.25)."""
    v = s.to_numpy(dtype=float, copy=True)
    bad = np.zeros(v.shape, bool)
    for nv in NULL_SENTINELS:
        bad |= np.isclose(v, nv, rtol=0, atol=1e-3) if abs(nv) < 1e6 else (np.abs(v) >= 1e29)
    v[bad] = np.nan
    return pd.Series(v, index=s.index, name=s.name), int(bad.sum())


def clip_physical(s: pd.Series, lo: float, hi: float) -> tuple[pd.Series, int]:
    v = s.to_numpy(dtype=float, copy=True)
    bad = (v < lo) | (v > hi)
    v[bad] = np.nan
    return pd.Series(v, index=s.index, name=s.name), int(np.nansum(bad))


def despike(s: pd.Series, window: int = 21, n_mad: float = 6.0, log_scale: bool = False,
            min_abs: float = 0.0) -> tuple[pd.Series, int]:
    """Null samples that sit more than ``n_mad`` robust standard deviations from a rolling median.

    ``min_abs`` is a floor on the deviation in curve units, so quiet intervals with a tiny MAD do not
    produce false spikes (e.g. 0.05 g/cm3 for RHOB). Resistivity should use ``log_scale=True``.
    """
    x = s.astype(float)
    if log_scale:
        x = np.log10(x.where(x > 0))
    med = x.rolling(window, center=True, min_periods=window // 3).median()
    resid = (x - med).abs()
    mad = resid.rolling(window * 5, center=True, min_periods=window).median() * 1.4826
    thresh = np.maximum(n_mad * mad, min_abs)
    spikes = (resid > thresh) & x.notna()
    out = s.astype(float).copy()
    out[spikes] = np.nan
    return out, int(spikes.sum())


def find_flatlines(s: pd.Series, min_length: float, tol: float = 1e-6) -> list[tuple[float, float, float]]:
    """Return (top, base, value) for runs where the curve is constant for at least ``min_length``."""
    v = s.dropna()
    if v.size < 3:
        return []
    same = np.abs(np.diff(v.to_numpy())) <= tol
    idx = v.index.to_numpy()
    out, start = [], None
    for i, flag in enumerate(np.append(same, False)):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            top, base = idx[start], idx[i]
            if base - top >= min_length:
                out.append((float(top), float(base), float(v.iloc[start])))
            start = None
    return out


def logged_below_td(well: Well) -> float:
    """Depth (in well units) by which any curve extends below the stated TD; 0 if none."""
    if well.logs is None or well.td is None or pd.isna(well.td):
        return 0.0
    valid = well.logs.dropna(how="all")
    if valid.empty:
        return 0.0
    return max(0.0, float(valid.index.max() - well.td))


def depth_match(ref: pd.Series, target: pd.Series, max_shift: float, step: float | None = None
                ) -> tuple[float, float]:
    """Estimate the bulk depth shift to apply to ``target`` so it lines up with ``ref``.

    Both curves are resampled to a common grid, standardised and cross-correlated over the overlap.
    Returns ``(shift, correlation)``; apply as ``target.index + shift``.
    """
    step = step or float(np.median(np.diff(ref.index)))
    top = max(ref.dropna().index.min(), target.dropna().index.min())
    base = min(ref.dropna().index.max(), target.dropna().index.max())
    if base - top < 10 * max_shift:
        raise ValueError("Overlap too short for depth matching")
    grid = np.arange(top, base, step)

    def prep(s):
        y = np.interp(grid, s.dropna().index, s.dropna().to_numpy())
        return (y - y.mean()) / (y.std() or 1.0)

    a, b = prep(ref), prep(target)
    n = int(round(max_shift / step))
    best, best_r = 0, -np.inf
    for k in range(-n, n + 1):
        if k >= 0:
            r = np.corrcoef(a[k:], b[: len(b) - k])[0, 1]
        else:
            r = np.corrcoef(a[: len(a) + k], b[-k:])[0, 1]
        if r > best_r:
            best, best_r = k, r
    return best * step, float(best_r)


def splice(runs: list[pd.DataFrame], priority: str = "shallow_first") -> pd.DataFrame:
    """Splice logging runs onto one depth grid.

    Where runs overlap the later run (``priority='deep_first'``) or the earlier run wins.
    Adds a ``RUN`` column recording which run each sample came from.
    """
    step = float(np.median(np.diff(runs[0].index)))
    top = min(r.index.min() for r in runs)
    base = max(r.index.max() for r in runs)
    grid = np.round(np.arange(top, base + step / 2, step), 4)
    cols = sorted(set().union(*[r.columns for r in runs]))
    out = pd.DataFrame(np.nan, index=grid, columns=cols)
    out["RUN"] = np.nan
    order = list(enumerate(runs))
    if priority == "shallow_first":
        order = order[::-1]
    for i, r in order:
        rr = r.reindex(grid, method="nearest", tolerance=step / 2)
        rr = rr.where(rr.index.to_series().between(r.index.min(), r.index.max()), axis=0)
        mask = rr.notna().any(axis=1)
        for c in r.columns:
            out.loc[rr[c].notna(), c] = rr.loc[rr[c].notna(), c]
        out.loc[mask, "RUN"] = i + 1
    out.index.name = runs[0].index.name
    return out


def condition_well(well: Well, lineage: Lineage, flat_min: float | None = None,
                   despike_cfg: dict | None = None) -> tuple[pd.DataFrame, list[Issue]]:
    """Run nulls -> physical limits -> despike -> flat-line checks on every curve of a well.

    Returns the conditioned logs and the issues found. The input well is not modified.
    """
    issues: list[Issue] = []
    out = well.logs.copy()
    step = abs(float(np.median(np.diff(out.index)))) if len(out) > 1 else 0.5
    flat_min = flat_min or (30.0 if well.depth_unit == "m" else 100.0)
    cfg = {"RHOB": dict(min_abs=0.08), "DT": dict(min_abs=8), "NPHI": dict(min_abs=0.06),
           "GR": dict(min_abs=25), "RDEEP": dict(log_scale=True, min_abs=0.4),
           "RMED": dict(log_scale=True, min_abs=0.4), "RSHAL": dict(log_scale=True, min_abs=0.4),
           "CALI": dict(min_abs=1.5), "DTS": dict(min_abs=15)}
    cfg.update(despike_cfg or {})
    for c in out.columns:
        canon = c.split("_ALT")[0]
        if canon not in PHYSICAL_RANGES:
            # unrecognised curves (bit size, flags, tool status, computed outputs) are passed through
            # untouched apart from null handling: conditioning rules only apply to known measurements
            out[c], _ = fix_nulls(out[c])
            continue
        s, n = fix_nulls(out[c])
        if n:
            sev = "major" if (s.dropna().size and n > 0) else "minor"
            issues.append(Issue(sev, "nulls", f"{n} undeclared null values (e.g. +999.25) set to null",
                                well.name, c))
            lineage.add("condition", "fix nulls", well.name, c, n)
        if canon in PHYSICAL_RANGES:
            lo, hi = PHYSICAL_RANGES[canon]
            s, n = clip_physical(s, lo, hi)
            if n:
                frac = n / max(1, out[c].notna().sum())
                issues.append(Issue("major" if frac > 0.01 else "minor", "range",
                                    f"{n} samples outside physical range [{lo}, {hi}] set to null",
                                    well.name, c))
                lineage.add("condition", "physical limits", well.name, c, n, lo=lo, hi=hi)
        if canon in cfg:
            s, n = despike(s, **cfg[canon])
            if n:
                issues.append(Issue("minor", "spike", f"{n} spikes removed", well.name, c))
                lineage.add("condition", "despike", well.name, c, n, **cfg[canon])
        for top, base, val in find_flatlines(s, flat_min):
            issues.append(Issue("major", "flatline",
                                f"Constant value {val:g} from {top:.1f} to {base:.1f} {well.depth_unit} "
                                f"({base - top:.0f} {well.depth_unit}); nulled as fill/filler", well.name, c))
            s[(s.index >= top) & (s.index <= base)] = np.nan
            lineage.add("condition", "null flat-line", well.name, c,
                        int(round((base - top) / step)), top=top, base=base, value=val)
        out[c] = s
    excess = logged_below_td(well)
    if excess > 2 * step and well.td is not None:
        issues.append(Issue("major", "depth", f"Logs extend {excess:.1f} {well.depth_unit} below stated TD "
                                              f"{well.td:g}; TD or depth reference inconsistent", well.name))
    return out, issues
