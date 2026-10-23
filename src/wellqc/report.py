"""Quality scoring and the Markdown audit report."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .model import Issue

SEVERITY_ORDER = {"critical": 0, "major": 1, "minor": 2, "info": 3}
KEY_CURVES = ("GR", "DT", "RHOB", "NPHI", "RDEEP", "DTS")


def curve_scores(raw: pd.DataFrame, clean: pd.DataFrame, well: str) -> pd.DataFrame:
    """Score each curve 0-100.

    * coverage (50): share of the well's logged interval where the conditioned curve has data
    * validity (30): share of the raw samples that survived nulls/limits/despike/flat-line checks
    * continuity (20): penalises fragmented curves (many short gaps inside the curve's own interval)
    """
    rows = []
    span = clean.dropna(how="all")
    n_well = max(1, len(span))
    for c in clean.columns:
        r, s = raw[c].dropna(), clean[c]
        if r.empty:
            rows.append({"well": well, "curve": c, "score": 0.0, "coverage": 0.0, "validity": 0.0,
                         "continuity": 0.0, "top": np.nan, "base": np.nan})
            continue
        cov = s.notna().sum() / n_well
        val = s.notna().sum() / max(1, r.size)
        inside = s.loc[r.index.min():r.index.max()]
        gaps = int((inside.notna().astype(int).diff() == -1).sum())
        cont = 1.0 / (1.0 + gaps / 10.0)
        rows.append({"well": well, "curve": c, "score": round(100 * (0.5 * cov + 0.3 * val + 0.2 * cont), 1),
                     "coverage": round(cov, 3), "validity": round(val, 3), "continuity": round(cont, 3),
                     "top": float(r.index.min()), "base": float(r.index.max())})
    return pd.DataFrame(rows)


def issues_frame(issues: list[Issue]) -> pd.DataFrame:
    df = pd.DataFrame([i.as_dict() for i in issues])
    if df.empty:
        return pd.DataFrame(columns=["severity", "category", "message", "well", "curve", "source"])
    df["_o"] = df["severity"].map(SEVERITY_ORDER)
    return df.sort_values(["_o", "category", "well"]).drop(columns="_o").reset_index(drop=True)


def _md_table(df: pd.DataFrame, max_rows: int = 60) -> str:
    if df is None or df.empty:
        return "_none_\n"
    df = df.head(max_rows).fillna("")
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(r[c]).replace("|", "/") for c in cols) + " |")
    return "\n".join(out) + "\n"


def write_report(out_dir: Path, inventory: pd.DataFrame, master: pd.DataFrame, duplicates: pd.DataFrame,
                 issues: pd.DataFrame, scores: pd.DataFrame, segy: pd.DataFrame, checkshots: pd.DataFrame,
                 prediction: pd.DataFrame | None, lineage: pd.DataFrame, figures: list[str]) -> Path:
    sev = issues["severity"].value_counts().reindex(["critical", "major", "minor", "info"]).fillna(0).astype(int)
    well_scores = (scores[scores["curve"].isin(KEY_CURVES)].groupby("well")["score"].mean().round(1)
                   .rename("key_curve_score").reset_index().sort_values("key_curve_score"))
    md = ["# Well data audit report\n",
          f"Files scanned: **{len(inventory)}** · physical wells: **{len(master)}** · "
          f"issues: **{sev['critical']} critical, {sev['major']} major, {sev['minor']} minor**\n",
          "## Well master (one row per physical well)\n", _md_table(master),
          "## Duplicate wells across projects\n", _md_table(duplicates),
          "## Issues\n", _md_table(issues[issues.severity != "info"], 200),
          "## Curve quality (key curves, mean score per well)\n", _md_table(well_scores),
          "## SEG-Y\n", _md_table(segy),
          "## Checkshots\n", _md_table(checkshots)]
    if prediction is not None and not prediction.empty:
        md += ["## Missing-curve prediction: blind-well validation\n", _md_table(prediction.round(3))]
    if figures:
        md.append("## Figures\n")
        md += [f"![{Path(f).stem}](figures/{Path(f).name})\n" for f in figures]
    md += ["## Scoring method\n",
           "Curve score = 50 x coverage of the well's logged interval + 30 x share of raw samples that passed "
           "null, physical-limit, despike and flat-line checks + 20 x continuity (1 / (1 + gaps/10)).\n",
           f"## Processing lineage\n{len(lineage)} recorded steps; full detail in `lineage.csv`.\n"]
    path = out_dir / "audit_report.md"
    path.write_text("\n".join(md))
    return path
