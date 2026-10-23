"""Checkshot / VSP quality control: interval velocities, outliers and quantisation-aware smoothing."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline

from .model import Issue


def interval_velocity(depth, twt_ms):
    """Interval velocity (m/s or ft/s, matching depth units) between consecutive levels."""
    d, t = np.asarray(depth, float), np.asarray(twt_ms, float)
    return 2.0 * np.diff(d) / (np.diff(t) / 1000.0)


def qc_checkshot(cs: pd.DataFrame, well: str | None = None, rel_tol: float = 0.15,
                 pick_sigma_ms: float = 1.0, time_quantum_ms: float | None = None
                 ) -> tuple[pd.DataFrame, list[Issue]]:
    """QC a checkshot table with columns DEPTH, TWT (ms).

    1. Levels whose time does not increase with depth are flagged ``non-monotonic``.
    2. An interval velocity deviates when it differs from the median of its neighbours by more than
       ``rel_tol`` *and* by more than TWT rounding plus normal pick noise (``pick_sigma_ms`` on each
       level) can explain over that interval. Two neighbouring intervals deviating in opposite
       directions identify a ``mis-picked level``; a lone deviating interval is an
       ``anomalous interval`` (thin layer or pick error) to verify.
    3. A smoothing spline is fitted to TWT(depth) through the unflagged levels, with the smoothing
       set by the expected pick/rounding error, giving a cleaned time-depth and interval-velocity
       profile (``TWT_SMOOTH``, ``VINT_SMOOTH``).
    """
    issues: list[Issue] = []
    df = cs.sort_values("DEPTH").reset_index(drop=True).copy()
    d, t = df["DEPTH"].to_numpy(float), df["TWT"].to_numpy(float)
    df["FLAG"] = ""
    nonmono = np.where(np.diff(t) <= 0)[0] + 1
    if nonmono.size:
        df.loc[nonmono, "FLAG"] = "non-monotonic"
        issues.append(Issue("critical", "checkshot", f"{nonmono.size} levels with non-increasing time", well=well))
    vint = pd.Series(np.r_[np.nan, interval_velocity(d, t)])
    df["VINT"] = vint
    if time_quantum_ms is None:
        time_quantum_ms = 1.0 if np.allclose(t, np.round(t)) else 0.0

    # neighbour median excluding the interval itself
    nb = pd.concat([vint.shift(k) for k in (-2, -1, 1, 2)], axis=1).median(axis=1)
    signed = vint / nb - 1
    rel = signed.abs()
    dt = pd.Series(np.r_[np.nan, np.diff(t)])
    # relative error a time difference can carry from rounding plus normal pick noise on both levels
    allowance = ((time_quantum_ms + 2 * pick_sigma_ms) / dt).fillna(0)
    big = ((rel > rel_tol) & (rel > allowance + 0.03)).to_numpy()
    quant = (rel > rel_tol) & ~pd.Series(big)
    sg = np.sign(signed.fillna(0).to_numpy())
    # VINT[i] spans levels i-1 -> i. A mis-picked level i makes VINT[i] and VINT[i+1] deviate in
    # opposite directions; a single deviating interval is more likely a thin layer (or a pick on
    # the last level) and is reported for review instead.
    mispick, anomalous = [], []
    i = 1
    while i < len(df):
        if big[i] and i + 1 < len(df) and big[i + 1] and sg[i] != sg[i + 1]:
            mispick.append(i)
            i += 2
            continue
        if big[i]:
            anomalous.append(i)
        i += 1
    for i in mispick:
        if df.at[i, "FLAG"] == "":
            df.at[i, "FLAG"] = "mis-picked level"
    for i in anomalous:
        if df.at[i, "FLAG"] == "":
            df.at[i, "FLAG"] = "anomalous interval"
    if quant.any():
        issues.append(Issue("minor", "checkshot",
                            f"{int(quant.sum())} interval velocities scatter > {rel_tol:.0%} but within what "
                            f"{time_quantum_ms:g} ms rounding and +/-{pick_sigma_ms:g} ms pick noise allow over short "
                            "intervals; use the smoothed profile", well=well))
    if mispick:
        lv = ", ".join(f"{d[i]:g}" for i in mispick[:6])
        issues.append(Issue("major", "checkshot", f"{len(mispick)} mis-picked level(s) at depth {lv}: the intervals "
                                                  "above and below deviate in opposite directions; excluded", well=well))
    if anomalous:
        lv = ", ".join(f"{d[i - 1]:g}-{d[i]:g}" for i in anomalous[:6])
        issues.append(Issue("minor", "checkshot", f"{len(anomalous)} anomalous interval velocit"
                                                  f"{'y' if len(anomalous) == 1 else 'ies'} ({lv}): thin layer or pick "
                                                  "error; verify against sonic", well=well))

    keep = (df["FLAG"] == "").to_numpy()
    sigma = max(pick_sigma_ms, time_quantum_ms / np.sqrt(12.0))
    spl = UnivariateSpline(d[keep], t[keep], k=3, s=keep.sum() * sigma ** 2)
    df["TWT_SMOOTH"] = spl(d)
    slope = spl.derivative()(d)  # ms per depth unit
    df["VINT_SMOOTH"] = 2000.0 / slope
    if (slope <= 0).any():
        issues.append(Issue("critical", "checkshot", "Smoothed time-depth curve is not monotonic", well=well))
    return df, issues
