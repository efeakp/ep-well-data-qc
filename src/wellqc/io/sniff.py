"""Identify E&P file types by content, not extension (legacy folders rarely use consistent names)."""

from __future__ import annotations

import re
import struct
from pathlib import Path


def sniff(path: str | Path) -> str:
    """Return one of: las, well_head, deviation, well_tops, checkshot, segy, surface, faults, unknown."""
    path = Path(path)
    if path.stat().st_size >= 3600 + 240:
        with open(path, "rb") as fh:
            head = fh.read(3600)
        fmt = struct.unpack(">h", head[3224:3226])[0]
        ns = struct.unpack(">h", head[3220:3222])[0]
        dt = struct.unpack(">h", head[3216:3218])[0]
        if fmt in (1, 2, 3, 5, 6, 8) and 0 < ns < 65000 and 0 < dt < 100000 and b"~V" not in head[:200]:
            return "segy"
    try:
        text = path.read_text(errors="replace")[:4000]
    except OSError:
        return "unknown"
    up = text.upper()
    first = up.lstrip().splitlines()[0] if up.strip() else ""
    if re.search(r"^~V", up, re.M):
        return "las"
    if "WELL TRACE FROM PETREL" in up:
        return "deviation"
    if "PETREL WELL TOPS" in up:
        return "well_tops"
    if "PETREL WELL HEAD" in up or (first.startswith("NAME") and "UWI" in first):
        return "well_head"
    if up.startswith("FSASCI"):
        return "surface"
    if up.startswith("FFASCI"):
        return "faults"
    lines = [ln for ln in up.splitlines() if ln.strip()][:6]
    if any(w in " ".join(lines[:2]) for w in ("TWT", "TIME", "CHECKSHOT", "VSP")):
        return "checkshot"
    return "unknown"
