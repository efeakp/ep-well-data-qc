"""Fast SEG-Y header scan: geometry, byte locations, truncation and dead traces.

Works on the raw bytes with ``numpy.memmap`` so it still reports on files that standard readers
refuse to open (e.g. truncated transfers or non-standard inline/crossline byte positions).
"""

from __future__ import annotations

import codecs
import re
import struct
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..model import Issue

_BYTES_PER_SAMPLE = {1: 4, 2: 4, 3: 2, 5: 4, 6: 8, 8: 1}
# Candidate 4-byte trace-header locations (1-based byte numbers) used for inline/crossline/coords
_CANDIDATES = {5: "TRACE_SEQUENCE_FILE", 9: "FIELD_RECORD", 17: "ENERGY_SOURCE_POINT", 21: "CDP",
               181: "CDP_X", 185: "CDP_Y", 189: "INLINE_3D", 193: "CROSSLINE_3D",
               73: "SOURCE_X", 77: "SOURCE_Y", 81: "GROUP_X", 85: "GROUP_Y"}


@dataclass
class SegyInfo:
    path: str
    text_header: str
    sample_format: int
    sample_interval_us: int
    n_samples: int
    start_time_ms: float
    n_traces_complete: int
    leftover_bytes: int
    inline_byte: int | None = None
    xline_byte: int | None = None
    x_byte: int | None = None
    y_byte: int | None = None
    inline_range: tuple[int, int] | None = None
    xline_range: tuple[int, int] | None = None
    n_inlines: int | None = None
    n_xlines: int | None = None
    expected_traces: int | None = None
    grid: dict | None = None  # affine X/Y = f(inline, xline)
    dead_traces: int = 0
    amplitude_pct: tuple[float, ...] = ()
    text_claims: dict = field(default_factory=dict)
    issues: list[Issue] = field(default_factory=list)

    @property
    def missing_traces(self) -> int:
        return max(0, (self.expected_traces or self.n_traces_complete) - self.n_traces_complete)


def _decode_text(raw: bytes) -> str:
    ebcdic = codecs.decode(raw, "cp500", errors="replace")
    ascii_ = raw.decode("latin1", errors="replace")
    # EBCDIC-decoded text of an ASCII header is mostly non-printable; pick the more readable one
    score = lambda s: sum(ch.isalnum() or ch == " " for ch in s)  # noqa: E731
    txt = ebcdic if score(ebcdic) > score(ascii_) else ascii_
    return "\n".join(txt[i:i + 80].rstrip() for i in range(0, 3200, 80))


def scan_segy(path: str | Path, n_amp_traces: int = 500) -> SegyInfo:
    path = Path(path)
    with open(path, "rb") as fh:
        head = fh.read(3600)
    text = _decode_text(head[:3200])
    bh = head[3200:3600]
    fmt = struct.unpack(">h", bh[24:26])[0]
    dt = struct.unpack(">h", bh[16:18])[0]
    ns = struct.unpack(">h", bh[20:22])[0]
    bps = _BYTES_PER_SAMPLE.get(fmt, 4)
    tl = 240 + ns * bps
    size = path.stat().st_size
    n = (size - 3600) // tl
    left = (size - 3600) - n * tl
    mm = np.memmap(path, dtype=np.uint8, mode="r", offset=3600, shape=(n * tl,)).reshape(n, tl)

    def i4(byte):
        o = byte - 1
        return mm[:, o:o + 4].copy().view(">i4").ravel()

    def i2(byte):
        o = byte - 1
        return mm[:, o:o + 2].copy().view(">i2").ravel()

    delay = i2(109)
    info = SegyInfo(str(path), text, fmt, dt, ns, float(np.median(delay)) if n else 0.0, n, left)
    iss = info.issues
    if fmt not in _BYTES_PER_SAMPLE:
        iss.append(Issue("critical", "segy", f"Unknown sample format code {fmt}", source=path.name))
    if left:
        iss.append(Issue("critical", "segy", f"File ends with a partial trace ({left} bytes); file is truncated",
                         source=path.name))

    # Text-header claims (Petrel writes First/Last inline and xline)
    claims = {}
    for key, pat in {"il_first": r"First inline:\s*(\d+)", "il_last": r"Last inline:\s*(\d+)",
                     "xl_first": r"First xline:\s*(\d+)", "xl_last": r"Last xline:\s*(\d+)"}.items():
        m = re.search(pat, text, re.I)
        if m:
            claims[key] = int(m.group(1))
    info.text_claims = claims

    # Detect which byte locations hold inline / crossline / coordinates
    vals = {b: i4(b) for b in _CANDIDATES}
    varying = {b: v for b, v in vals.items() if len(np.unique(v)) > 1}
    small = {b: v for b, v in varying.items() if np.abs(v).max() < 1_000_000 and b not in (73, 77, 81, 85, 181, 185)}
    big = {b: v for b, v in varying.items() if b in (73, 77, 81, 85, 181, 185)}
    if len(small) >= 2:
        # inline changes slowest along the file; crossline fastest
        changes = {b: int((np.diff(v) != 0).sum()) for b, v in small.items()}
        info.inline_byte = min(changes, key=changes.get)
        # crossline: changes often but repeats (unlike a running trace counter, which is all-unique)
        rest = {b: len(np.unique(v)) for b, v in small.items() if b != info.inline_byte}
        info.xline_byte = min(rest, key=rest.get)
    if 189 in small and 193 in small:
        info.inline_byte, info.xline_byte = 189, 193
    for xb, yb in ((181, 185), (73, 77), (81, 85)):
        if xb in big and yb in big:
            info.x_byte, info.y_byte = xb, yb
            break
    if info.inline_byte and (info.inline_byte, info.xline_byte) != (189, 193):
        iss.append(Issue("major", "segy",
                         f"Inline/crossline stored at non-standard bytes {info.inline_byte}/{info.xline_byte} "
                         f"({_CANDIDATES[info.inline_byte]}/{_CANDIDATES[info.xline_byte]}); bytes 189/193 are "
                         "empty. Set these byte locations explicitly on load.", source=path.name))
    if info.x_byte and info.x_byte != 181:
        iss.append(Issue("minor", "segy", f"Coordinates stored in bytes {info.x_byte}/{info.y_byte} "
                                          f"({_CANDIDATES[info.x_byte]}), not CDP_X/CDP_Y (181/185)",
                         source=path.name))

    if info.inline_byte:
        il, xl = vals[info.inline_byte], vals[info.xline_byte]
        il_first, il_last = claims.get("il_first", il.min()), claims.get("il_last", il.max())
        xl_first, xl_last = claims.get("xl_first", xl.min()), claims.get("xl_last", xl.max())
        info.inline_range, info.xline_range = (int(il_first), int(il_last)), (int(xl_first), int(xl_last))
        info.n_inlines = int(il_last - il_first + 1)
        info.n_xlines = int(xl_last - xl_first + 1)
        info.expected_traces = info.n_inlines * info.n_xlines
        if info.missing_traces:
            iss.append(Issue("critical", "segy",
                             f"{info.missing_traces} of {info.expected_traces} traces missing "
                             f"(last complete trace IL {il[-1]} XL {xl[-1]}); re-export from source",
                             source=path.name))
        if info.x_byte:
            scal = i2(71)
            s = np.where(scal < 0, -1.0 / scal, np.where(scal > 0, scal, 1.0)).astype(float)
            x, y = vals[info.x_byte] * s, vals[info.y_byte] * s
            A = np.column_stack([np.ones(n), il, xl])
            cx, *_ = np.linalg.lstsq(A, x, rcond=None)
            cy, *_ = np.linalg.lstsq(A, y, rcond=None)
            info.grid = {"x0": cx[0], "dx_il": cx[1], "dx_xl": cx[2],
                         "y0": cy[0], "dy_il": cy[1], "dy_xl": cy[2],
                         "il_spacing": float(np.hypot(cx[1], cy[1])), "xl_spacing": float(np.hypot(cx[2], cy[2])),
                         "xmin": float(x.min()), "xmax": float(x.max()), "ymin": float(y.min()), "ymax": float(y.max())}

    # Amplitude statistics and dead traces
    data = mm[:, 240:]
    if fmt == 5:
        samp = data.copy().view(">f4").reshape(n, ns)
        info.dead_traces = int((np.abs(samp).max(axis=1) == 0).sum())
        pick = samp[np.linspace(0, n - 1, min(n, n_amp_traces)).astype(int)]
        info.amplitude_pct = tuple(np.round(np.percentile(pick, [0, 1, 50, 99, 100]), 3))
        if info.dead_traces:
            iss.append(Issue("minor", "segy", f"{info.dead_traces} dead (all-zero) traces", source=path.name))
    return info


def grid_to_xy(grid: dict, il, xl):
    il, xl = np.asarray(il, float), np.asarray(xl, float)
    return (grid["x0"] + grid["dx_il"] * il + grid["dx_xl"] * xl,
            grid["y0"] + grid["dy_il"] * il + grid["dy_xl"] * xl)


def xy_to_grid(grid: dict, x, y):
    M = np.array([[grid["dx_il"], grid["dx_xl"]], [grid["dy_il"], grid["dy_xl"]]])
    rhs = np.vstack([np.asarray(x, float) - grid["x0"], np.asarray(y, float) - grid["y0"]])
    il, xl = np.linalg.solve(M, rhs)
    return il, xl


def same_grid(a: SegyInfo, b: SegyInfo, tol_m: float = 1.0) -> bool:
    """True when two surveys share the same inline/crossline-to-XY transform (e.g. crops of one cube)."""
    if not (a.grid and b.grid):
        return False
    il = np.array([a.inline_range[0], a.inline_range[1]])
    xl = np.array([a.xline_range[0], a.xline_range[1]])
    xa, ya = grid_to_xy(a.grid, il, xl)
    xb, yb = grid_to_xy(b.grid, il, xl)
    return bool(np.all(np.hypot(xa - xb, ya - yb) < tol_m))
