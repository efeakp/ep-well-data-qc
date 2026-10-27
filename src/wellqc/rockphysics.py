"""Rock physics: empirical velocity/density relations, Gassmann fluid substitution, a simple
rock physics template (RPT), and machine-learning prediction of missing elastic logs.

Units: velocities in km/s, densities in g/cm3, moduli in GPa unless stated otherwise.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

US_PER_FT_TO_KMS = 304.8  # Vp[km/s] = 304.8 / DT[us/ft]


def dt_to_vp(dt_usft):
    return US_PER_FT_TO_KMS / np.asarray(dt_usft, dtype=float)


def vp_to_dt(vp_kms):
    return US_PER_FT_TO_KMS / np.asarray(vp_kms, dtype=float)


# --------------------------------------------------------------------------- empirical relations
def gardner_rho(vp_kms, a: float = 1.741, b: float = 0.25):
    """Gardner et al. (1974): rho = a * Vp^b, with Vp in km/s (a = 0.31 for m/s)."""
    return a * np.asarray(vp_kms, float) ** b


def gardner_vp(rho, a: float = 1.741, b: float = 0.25):
    return (np.asarray(rho, float) / a) ** (1 / b)


def castagna_mudrock_vs(vp_kms):
    """Castagna et al. (1985) mudrock line: Vs = 0.8621 Vp - 1.1724 (km/s)."""
    return 0.8621 * np.asarray(vp_kms, float) - 1.1724


# Greenberg & Castagna (1992) linear Vp-Vs coefficients for brine-saturated end members (km/s)
GC_COEFFS = {"sandstone": (0.0, 0.80416, -0.85588), "limestone": (-0.05508, 1.01677, -1.03049),
             "dolomite": (0.0, 0.58321, -0.07775), "shale": (0.0, 0.76969, -0.86735)}


def greenberg_castagna_vs(vp_kms, fractions: dict[str, np.ndarray | float]):
    """Greenberg-Castagna Vs for a lithology mixture (brine saturated).

    ``fractions`` maps lithology -> volume fraction (e.g. {'sandstone': 1 - vsh, 'shale': vsh}).
    Uses the average of the Voigt (arithmetic) and Reuss (harmonic) bounds of the end-member Vs.
    """
    vp = np.asarray(vp_kms, float)
    tot = sum(np.asarray(f, float) for f in fractions.values())
    voigt, reuss = 0.0, 0.0
    for lith, f in fractions.items():
        a2, a1, a0 = GC_COEFFS[lith]
        vs_i = a2 * vp ** 2 + a1 * vp + a0
        w = np.asarray(f, float) / tot
        voigt = voigt + w * vs_i
        reuss = reuss + w / vs_i
    return 0.5 * (voigt + 1.0 / reuss)


# --------------------------------------------------------------------------- Gassmann
def moduli(vp, vs, rho):
    """Return (K, mu) in GPa from km/s and g/cm3."""
    mu = rho * vs ** 2
    return rho * vp ** 2 - 4.0 / 3.0 * mu, mu


def velocities(k, mu, rho):
    return np.sqrt((k + 4.0 / 3.0 * mu) / rho), np.sqrt(mu / rho)


def gassmann_kdry(k_sat, k_min, k_fl, phi):
    num = k_sat * (phi * k_min / k_fl + 1 - phi) - k_min
    den = phi * k_min / k_fl + k_sat / k_min - 1 - phi
    return num / den


def gassmann_ksat(k_dry, k_min, k_fl, phi):
    return k_dry + (1 - k_dry / k_min) ** 2 / (phi / k_fl + (1 - phi) / k_min - k_dry / k_min ** 2)


def fluid_substitute(vp, vs, rho, phi, k_fl1, rho_fl1, k_fl2, rho_fl2, k_min=37.0):
    """Gassmann substitution from fluid 1 to fluid 2. Returns (vp2, vs2, rho2)."""
    k1, mu = moduli(vp, vs, rho)
    kd = gassmann_kdry(k1, k_min, k_fl1, phi)
    k2 = gassmann_ksat(kd, k_min, k_fl2, phi)
    rho2 = rho + phi * (rho_fl2 - rho_fl1)
    vp2, vs2 = velocities(k2, mu, rho2)
    return vp2, vs2, rho2


# --------------------------------------------------------------------------- template
def rpt_lines(phi=None, vp_min: float = 5.95, rho_min: float = 2.65, vp_fl: float = 1.6,
              rho_brine: float = 1.05, k_brine: float = 2.8, rho_gas: float = 0.2, k_gas: float = 0.1,
              shale_vp=None) -> dict[str, pd.DataFrame]:
    """Build simple template trends in AI vs Vp/Vs space.

    * brine sand: Raymer-Hunt-Gardner Vp with Greenberg-Castagna sandstone Vs
    * gas sand: Gassmann substitution of the brine trend
    * shale: Castagna mudrock line with Gardner density
    Each frame has columns phi (or vp), vp, vs, rho, AI, VpVs.
    """
    phi = np.linspace(0.05, 0.35, 31) if phi is None else np.asarray(phi, float)
    vp_b = (1 - phi) ** 2 * vp_min + phi * vp_fl
    vs_b = greenberg_castagna_vs(vp_b, {"sandstone": 1.0})
    rho_b = (1 - phi) * rho_min + phi * rho_brine
    vp_g, vs_g, rho_g = fluid_substitute(vp_b, vs_b, rho_b, phi, k_brine, rho_brine, k_gas, rho_gas)
    shale_vp = np.linspace(2.0, 4.0, 21) if shale_vp is None else np.asarray(shale_vp, float)
    vs_s = castagna_mudrock_vs(shale_vp)
    rho_s = gardner_rho(shale_vp)

    def frame(**kw):
        df = pd.DataFrame(kw)
        df["AI"] = df.vp * df.rho
        df["VpVs"] = df.vp / df.vs
        return df

    return {"brine sand": frame(phi=phi, vp=vp_b, vs=vs_b, rho=rho_b),
            "gas sand": frame(phi=phi, vp=vp_g, vs=vs_g, rho=rho_g),
            "shale": frame(vp=shale_vp, vs=vs_s, rho=rho_s)}


def rpt_distance(ai, vpvs, template: dict[str, pd.DataFrame], tol: float = 0.25) -> np.ndarray:
    """Fraction-normalised distance from each sample to the nearest template trend (0 = on a trend).

    Samples with distance > ``tol`` are geologically implausible for the template and are flagged.
    """
    pts = np.column_stack([np.asarray(ai, float) / 10.0, np.asarray(vpvs, float) / 2.0])
    best = np.full(len(pts), np.inf)
    for df in template.values():
        line = np.column_stack([df.AI / 10.0, df.VpVs / 2.0])
        d = np.min(np.linalg.norm(pts[:, None, :] - line[None, :, :], axis=2), axis=1)
        best = np.minimum(best, d)
    return best


# --------------------------------------------------------------------------- ML prediction
def predict_curve(wells: dict[str, pd.DataFrame], target: str, features: list[str],
                  min_samples: int = 500, random_state: int = 0):
    """Predict ``target`` from ``features`` with leave-one-well-out (blind well) validation.

    ``wells`` maps well name -> logs frame (canonical mnemonics, plus optional 'TVDSS').
    Wells with the target are used for training/validation; the final model is trained on all of
    them and applied to every well. Returns (predictions dict, validation frame, model,
    blind predictions dict) where the blind predictions come from models that never saw that well.
    """
    from sklearn.ensemble import HistGradientBoostingRegressor

    def xy(df):
        X = df.reindex(columns=features)
        y = df[target] if target in df else pd.Series(np.nan, index=df.index)
        return X, y

    train = {n: d for n, d in wells.items() if target in d and d[target].notna().sum() >= min_samples}
    if len(train) < 2:
        raise ValueError(f"Need at least two wells with {target} for blind-well validation")
    # keep only features that every training well actually has where the target is measured;
    # otherwise a held-out well can leave the model with an all-empty input column
    requested = list(features)
    features = [f for f in requested
                if all(f in d and (d[target].notna() & d[f].notna()).sum() >= min_samples for d in train.values())]
    if len(features) < 2:
        raise ValueError(f"Fewer than two of {requested} are logged alongside {target} in every training well")
    rows, blind_preds = [], {}
    for blind in train:
        Xs, ys = zip(*(xy(d) for n, d in train.items() if n != blind), strict=True)
        X, y = pd.concat(Xs), pd.concat(ys)
        ok = y.notna() & X.notna().sum(axis=1).ge(len(features) - 1)
        m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, random_state=random_state)
        m.fit(X[ok], y[ok])
        Xb, yb = xy(train[blind])
        okb = yb.notna() & Xb.notna().sum(axis=1).ge(len(features) - 1)
        pred = m.predict(Xb[okb])
        blind_preds[blind] = pd.Series(pred, index=yb[okb].index, name=f"{target}_BLIND")
        err = pred - yb[okb].to_numpy()
        rows.append({"blind_well": blind, "features": "+".join(features), "n": int(okb.sum()), "rmse": float(np.sqrt(np.mean(err ** 2))),
                     "mae": float(np.mean(np.abs(err))),
                     "r2": float(1 - np.sum(err ** 2) / np.sum((yb[okb] - yb[okb].mean()) ** 2))})
    X, y = zip(*(xy(d) for d in train.values()), strict=True)
    X, y = pd.concat(X), pd.concat(y)
    ok = y.notna() & X.notna().sum(axis=1).ge(len(features) - 1)
    model = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, random_state=random_state)
    model.fit(X[ok], y[ok])
    preds = {}
    for n, d in wells.items():
        Xd, _ = xy(d)
        okd = Xd.notna().sum(axis=1).ge(len(features) - 1)
        p = pd.Series(np.nan, index=d.index, name=f"{target}_PRED")
        if okd.any():
            p[okd] = model.predict(Xd[okd])
        preds[n] = p
    return preds, pd.DataFrame(rows), model, blind_preds
