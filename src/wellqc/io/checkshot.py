"""Checkshot / VSP time-depth file reader that copes with unlabelled or swapped columns."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from ..model import Issue

_TIME_WORDS = ("twt", "time", "ms", "owt")
_DEPTH_WORDS = ("depth", "md", "tvd", "z")


def read_checkshot(path: str | Path) -> tuple[dict, pd.DataFrame, list[Issue]]:
    """Read a two-column time-depth table.

    Column roles are decided from header words (``TWT``, ``2-way time``, ``Measured Depth``...).
    When the header is missing or ambiguous the file is read as DEPTH, TWT and flagged for review
    rather than guessed silently.
    Returns ``(meta, frame[DEPTH, TWT], issues)`` with TWT in ms.
    """
    path = Path(path)
    lines = [ln.rstrip() for ln in path.read_text(errors="replace").splitlines() if ln.strip()]
    issues: list[Issue] = []
    meta: dict = {"source": path.name, "well": None, "depth_ref": None}
    header_tokens: list[str] = []
    rows = []
    for ln in lines:
        tok = re.split(r"\t|\s{2,}|,", ln.strip())
        nums = []
        for t in tok:
            try:
                nums.append(float(t))
            except ValueError:
                nums = None
                break
        if nums and len(nums) >= 2:
            rows.append(nums[:2])
        elif not rows:
            header_tokens.append(ln.strip())
    if not rows:
        raise ValueError(f"{path.name}: no numeric time-depth rows")
    a = np.array(rows)

    # Header analysis: last header line holds column names; an earlier one may hold the well name
    col_line = header_tokens[-1] if header_tokens else ""
    cols = [c.strip().lower() for c in re.split(r"\t|\s{2,}|,", col_line) if c.strip()]
    if len(header_tokens) > 1:
        meta["well"] = re.split(r"\t|\s{2,}", header_tokens[0])[0].strip()
    elif cols and not any(w in cols[0] for w in _TIME_WORDS + _DEPTH_WORDS):
        meta["well"] = cols[0]

    def role(c: str) -> str | None:
        if any(w in c for w in ("twt", "2-way", "two-way", "time")):
            return "time"
        if any(w in c for w in ("depth", "md", "tvd")):
            return "depth"
        return None

    roles = [role(c) for c in cols[:2]] if len(cols) >= 2 else [None, None]
    if roles == ["time", "depth"]:
        a = a[:, ::-1]
        issues.append(Issue("info", "format", "Columns ordered TIME, DEPTH; swapped to DEPTH, TWT",
                            source=path.name))
    elif roles != ["depth", "time"]:
        issues.append(Issue("major", "format", f"Could not identify time/depth columns from header {cols}; "
                                               "assumed DEPTH, TWT", source=path.name))
    joined = " ".join(cols)
    if "md" in joined or "measured" in joined:
        meta["depth_ref"] = "MD"
    elif "tvdss" in joined:
        meta["depth_ref"] = "TVDSS"
    elif "tvd" in joined:
        meta["depth_ref"] = "TVD"
    else:
        issues.append(Issue("major", "datum", "Checkshot depth reference (MD/TVD/TVDSS) and datum not stated",
                            well=meta["well"], source=path.name))
    df = pd.DataFrame(a, columns=["DEPTH", "TWT"])
    return meta, df, issues
