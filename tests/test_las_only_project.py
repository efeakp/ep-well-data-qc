"""A project made only of LAS files: no well header, bottom-up recording, '.F' depth units,
empty WELL field, metadata only in the LAS header, and reference curves (bit size) that must
not be treated as filler."""

import lasio
import numpy as np
import pytest

from wellqc.pipeline import run_audit


def _las(path, well, depth, curves, unit="M", params=None, well_items=None):
    las = lasio.LASFile()
    las.well["WELL"].value = well
    for k, (v, u) in (well_items or {}).items():
        las.well.append(lasio.HeaderItem(k, unit=u, value=v))
    for k, (v, u) in (params or {}).items():
        las.params.append(lasio.HeaderItem(k, unit=u, value=v))
    las.append_curve("DEPT", depth, unit=unit)
    for name, (vals, u) in curves.items():
        las.append_curve(name, vals, unit=u)
    las.write(str(path), version=2.0)


@pytest.fixture(scope="module")
def las_only(tmp_path_factory):
    root = tmp_path_factory.mktemp("las_only")
    rng = np.random.default_rng(3)
    z = np.arange(500.0, 1500.0, 0.1)
    gr = 60 + 40 * np.sin(z / 37) + rng.normal(0, 3, z.size)
    rhob = 2.2 + 0.0002 * (z - 500) + rng.normal(0, 0.01, z.size)
    dt = 140 - 0.04 * (z - 500) - 0.15 * (gr - 60) + rng.normal(0, 1.5, z.size)
    # well A: bottom-up, metadata in ~Parameter, constant bit-size curve
    zr = z[::-1]
    _las(root / "A-01.las", "ALPHA-01", zr, {"GR": (gr[::-1], "GAPI"), "RHOB": (rhob[::-1], "G/C3"),
                                             "DT": (dt[::-1], "US/F"), "BS": (np.full(z.size, 8.5), "IN")},
         params={"TDD": (1500.5, "M"), "APD": (7.9, "M")})
    # well B: in feet ('.F'), empty WELL field, X/Y in ~Well
    zf = z / 0.3048
    _las(root / "B-02.las", "", zf, {"GR": (gr + 5, "GAPI"), "RHOB": (rhob, "G/C3"), "DT": (dt + 2, "US/F")},
         unit="F", well_items={"X": (141000.0, "m"), "Y": (456000.0, "m"), "EKB": (30.0, "F")})
    # well C: normal, no DT (prediction target)
    _las(root / "C-03.las", "CHARLIE-03", z, {"GR": (gr - 5, "GAPI"), "RHOB": (rhob + 0.02, "G/C3")})
    return run_audit({"LAS_ONLY": root}, tmp_path_factory.mktemp("out"), figures=False)


def test_wells_built_from_las_headers(las_only):
    m = las_only.master.set_index("standard_name")
    assert set(m.index) == {"ALPHA-01", "B-02", "CHARLIE-03"}
    assert m.loc["ALPHA-01", "td"] == pytest.approx(1500.5)
    assert m.loc["ALPHA-01", "kb"] == pytest.approx(7.9)
    assert m.loc["B-02", "x"] == pytest.approx(141000.0)
    assert m.loc["B-02", "kb"] == pytest.approx(30 * 0.3048, abs=0.01)  # EKB given in feet


def test_feet_converted_and_bottom_up_sorted(las_only):
    b = las_only.clean_logs["B-02"]
    assert b.index.min() == pytest.approx(500.0, abs=0.01) and b.index.max() < 1501
    a = las_only.clean_logs["ALPHA-01"]
    assert a.index.is_monotonic_increasing
    msgs = [i.message for i in las_only.issues]
    assert any("bottom-up" in m for m in msgs)
    assert any("LAS WELL field is empty" in m for m in msgs)
    assert not any(i.severity == "critical" and i.category == "units" for i in las_only.issues)


def test_reference_curves_not_treated_as_filler(las_only):
    assert not any(i.category == "flatline" and i.curve == "BS" for i in las_only.issues)
    assert las_only.clean_logs["ALPHA-01"]["BS"].notna().all()


def test_prediction_uses_shared_features(las_only):
    v = las_only.validation
    assert not v.empty and set(v.features) == {"GR+RHOB+TVDSS"}
    assert las_only.predictions["CHARLIE-03"]["DT"].notna().sum() > 5000
