"""End-to-end: every defect injected by wellqc.synthetic must be found by the audit."""

from conftest import messages


def test_physical_well_count_and_duplicate(audit):
    assert len(audit.wells) == 7
    assert len(audit.master) == 6
    dup = audit.duplicates[audit.duplicates.verdict == "duplicate"]
    assert {frozenset((r.well_a, r.well_b)) for _, r in dup.iterrows()} == {frozenset(("ALPHA-03", "BETA 7"))}
    row = audit.master[audit.master.aliases.str.contains("BETA 7")].iloc[0]
    assert "ALPHA-03" in row.aliases and row.n_sources == 2


def test_sign_flipped_nulls_and_spikes(audit):
    assert any("+999.25" in m for m in messages(audit, "nulls", "ALPHA-02"))
    clean = audit.clean_logs["ALPHA-02"]
    assert clean["DT"].max() < 200
    assert clean["RHOB"].dropna().between(1.9, 2.9).all()


def test_flatline_filler(audit):
    msgs = messages(audit, "flatline", "ALPHA-01")
    assert msgs and "77.94" in msgs[0]


def test_logs_below_td(audit):
    assert any("below stated TD" in m for m in messages(audit, "depth", "ALPHA-02"))


def test_name_matching(audit):
    naming = " ".join(messages(audit, "naming"))
    assert "'ALPHAA-02'" in naming and "similar name" in naming
    assert "'CORE-04'" in naming and "wellhead location" in naming


def test_units(audit):
    units = messages(audit, "units")
    assert sum("NPHI labelled 'm3/m3' but values are in percent" in m for m in units) == 3
    assert any("VSH labelled '%' but values are fractions" in m for m in units)
    assert audit.clean_logs["BETA 5"]["NPHI"].max() < 1.0


def test_depth_match_and_splice(audit):
    msg = [m for m in messages(audit, "depth", "BETA 6") if "off depth" in m]
    assert msg, "multi-run depth shift not detected"
    shift = float(msg[0].split("off depth by ")[1].split()[0])
    assert abs(shift - 2.4) < 0.2


def test_missing_sonic_predicted(audit):
    pred = audit.predictions["ALPHA-04"]["DT"]
    assert pred.notna().sum() > 10_000
    assert audit.validation.r2.min() > 0.8


def test_checkshot(audit):
    cs = " ".join(messages(audit, "checkshot", "BETA 5"))
    assert "mis-picked level(s) at depth 1700" in cs
    assert any("swapped" in m for m in messages(audit, "format", "BETA 5"))


def test_segy(audit):
    segy = " ".join(messages(audit, "segy"))
    assert "non-standard bytes 5/21" in segy
    assert "37 of 2576 traces missing" in segy
    assert "partial trace" in segy


def test_outputs_written(audit):
    out = audit.report_path.parent
    for f in ("inventory.csv", "well_master.csv", "issues.csv", "curve_scores.csv", "lineage.csv",
              "audit_report.md"):
        assert (out / f).exists()
    assert len(list((out / "las").glob("*.las"))) == 6
