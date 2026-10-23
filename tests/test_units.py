"""Unit tests for the building blocks."""

import numpy as np
import pandas as pd
import pytest

from wellqc.checkshots import qc_checkshot
from wellqc.condition import depth_match, despike, find_flatlines, fix_nulls, splice
from wellqc.rockphysics import (
    castagna_mudrock_vs,
    fluid_substitute,
    gardner_rho,
    greenberg_castagna_vs,
    moduli,
    rpt_lines,
)
from wellqc.standardise import canonical_name, match_name, name_key


@pytest.mark.parametrize("name,key", [("DELTA-02", ("DELTA", 2, "")), ("BETA 05", ("BETA", 5, "")),
                                      ("omega_11", ("OMEGA", 11, "")), ("ABC-3ST1", ("ABC", 3, "ST1"))])
def test_name_key(name, key):
    assert name_key(name) == key


def test_match_name_fuzzy_requires_same_number():
    assert match_name("ALPHAA 2", ["ALPHA 2", "ALPHA 12"]) == ("ALPHA 2", "fuzzy")
    assert match_name("BETA 5", ["BETA 05"]) == ("BETA 05", "exact")
    assert match_name("CORE-04", ["DELTA-04"])[0] is None
    assert canonical_name("beta 5") == "BETA-05"


def test_fix_nulls_catches_sign_flipped_sentinel():
    s = pd.Series([1.0, 999.25, -999.25, 2.0, -9999.0])
    out, n = fix_nulls(s)
    assert n == 3 and out.isna().sum() == 3


def test_despike_removes_spikes_keeps_steps():
    rng = np.random.default_rng(0)
    x = np.r_[np.full(300, 2.3), np.full(300, 2.55)] + rng.normal(0, 0.01, 600)
    x[[50, 400]] += [0.8, -0.9]
    out, n = despike(pd.Series(x), min_abs=0.08)
    assert n == 2
    assert out.iloc[290:310].notna().sum() >= 19  # the real step is preserved


def test_flatlines():
    s = pd.Series(np.r_[np.arange(10.0), np.full(100, 5.0), np.arange(10.0)], index=np.arange(120) * 0.5)
    runs = find_flatlines(s, min_length=30)
    assert len(runs) == 1 and runs[0][2] == 5.0


def test_depth_match_recovers_shift():
    rng = np.random.default_rng(1)
    z = np.arange(1000, 1600, 0.1524)
    gr = pd.Series(np.convolve(rng.normal(0, 1, z.size), np.ones(25) / 25, "same"), index=z)
    shifted = pd.Series(gr.to_numpy(), index=z + 1.8)
    shift, r = depth_match(gr, shifted, max_shift=5)
    assert abs(shift + 1.8) < 0.16 and r > 0.95


def test_splice_prefers_shallow_run():
    a = pd.DataFrame({"GR": 1.0}, index=np.round(np.arange(0, 10, 0.5), 4))
    b = pd.DataFrame({"GR": 2.0}, index=np.round(np.arange(8, 20, 0.5), 4))
    out = splice([a, b])
    assert out.loc[9.0, "GR"] == 1.0 and out.loc[15.0, "GR"] == 2.0 and out["RUN"].iloc[-1] == 2


def test_rock_physics_relations():
    assert gardner_rho(3.0) == pytest.approx(2.29, abs=0.01)
    assert castagna_mudrock_vs(3.0) == pytest.approx(1.4139, abs=1e-4)
    assert greenberg_castagna_vs(3.0, {"sandstone": 1.0}) == pytest.approx(0.80416 * 3 - 0.85588)


def test_gassmann_round_trip():
    vp, vs, rho, phi = 3.2, 1.75, 2.25, 0.25
    vp2, vs2, rho2 = fluid_substitute(vp, vs, rho, phi, 2.8, 1.05, 0.1, 0.2)
    vp3, vs3, rho3 = fluid_substitute(vp2, vs2, rho2, phi, 0.1, 0.2, 2.8, 1.05)
    assert vp2 < vp and vs2 > vs  # gas lowers Vp, raises Vs (lower density)
    assert (vp3, vs3, rho3) == pytest.approx((vp, vs, rho), rel=1e-9)
    assert moduli(vp, vs, rho)[1] == pytest.approx(rho * vs ** 2)


def test_rpt_gas_sand_has_lower_vpvs():
    t = rpt_lines()
    assert (t["gas sand"].VpVs < t["brine sand"].VpVs).all()


def test_checkshot_mispick_vs_thin_layer():
    depth = np.arange(500.0, 3000.0, 50.0)
    v = np.full(depth.size, 3000.0)
    v[10:12] = 3600.0  # a real fast layer two intervals thick
    twt = np.round(np.r_[0, np.cumsum(2 * np.diff(depth) / v[1:] * 1000)] + 500)
    twt[30] += 8  # one mis-picked level
    df, issues = qc_checkshot(pd.DataFrame({"DEPTH": depth, "TWT": twt}))
    assert df.loc[30, "FLAG"] == "mis-picked level"
    assert (df.FLAG == "mis-picked level").sum() == 1
    assert df["VINT_SMOOTH"].between(2700, 3800).all()
