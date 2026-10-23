"""Readers for Petrel ASCII exports: well heads, deviation traces and well tops.

Petrel exports are loosely specified, so each reader returns its data together with a list of
:class:`~wellqc.model.Issue` objects describing anything inconsistent it noticed while parsing.
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path

import numpy as np
import pandas as pd

from ..model import Issue

_UNIT_RE = re.compile(r"#\s*Unit in (X and Y direction|depth):\s*(\w+)", re.I)


def _numeric(v: str):
    try:
        f = float(v)
    except ValueError:
        return v
    return np.nan if f == -999 or f == -999.25 else f


def _parse_header_block(lines: list[str]) -> tuple[list[str], int]:
    """Return the field names between BEGIN HEADER / END HEADER and the index after END HEADER."""
    start = next(i for i, ln in enumerate(lines) if ln.strip().upper() == "BEGIN HEADER")
    end = next(i for i, ln in enumerate(lines) if ln.strip().upper() == "END HEADER")
    return [ln.strip() for ln in lines[start + 1:end]], end + 1


def read_well_head(path: str | Path) -> tuple[pd.DataFrame, list[Issue]]:
    """Read a Petrel well-head export or a simple tab-delimited well header table.

    Returns a frame with at least ``name, uwi, x, y, kb, td`` plus ``xy_unit`` and ``depth_unit``.
    """
    path = Path(path)
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    issues: list[Issue] = []
    units = {"xy_unit": None, "depth_unit": None}
    for m in _UNIT_RE.finditer(text):
        units["xy_unit" if m.group(1).lower().startswith("x") else "depth_unit"] = m.group(2).lower()

    if any(ln.strip().upper() == "BEGIN HEADER" for ln in lines):
        fields, i0 = _parse_header_block(lines)
        rows = []
        for ln in lines[i0:]:
            if not ln.strip():
                continue
            tok = shlex.split(ln)
            if len(tok) != len(fields):
                issues.append(Issue("minor", "format",
                                    f"Row has {len(tok)} values but header defines {len(fields)} fields",
                                    well=tok[0] if tok else None, source=path.name))
            rows.append(tok[: len(fields)] + [""] * max(0, len(fields) - len(tok)))
        raw = pd.DataFrame(rows, columns=_dedupe(fields))
        col = {c.lower(): c for c in raw.columns}
        df = pd.DataFrame({
            "name": raw[col["name"]],
            "uwi": raw[col["uwi"]] if "uwi" in col else None,
            "x": raw[col["surface x"]].map(_numeric),
            "y": raw[col["surface y"]].map(_numeric),
            "kb": raw[col["kb"]].map(_numeric) if "kb" in col else np.nan,
            "td": raw[col["td (md)"]].map(_numeric) if "td (md)" in col else np.nan,
        })
    else:
        # Simple delimited table: Name UWI X Y KB TD (header names are often ambiguous)
        raw = pd.read_csv(path, sep=r"\t|\s{2,}|,", engine="python")
        cols = list(raw.columns)
        if len(cols) < 6:
            raise ValueError(f"{path.name}: cannot interpret well header with columns {cols}")
        if len(set(c.split(".")[0] for c in cols[2:4])) == 1:
            issues.append(Issue("minor", "format",
                                f"Ambiguous coordinate column names {cols[2:4]}; assumed X then Y",
                                source=path.name))
        df = pd.DataFrame({"name": raw.iloc[:, 0].astype(str).str.strip(),
                           "uwi": raw.iloc[:, 1].astype(str).str.strip(),
                           "x": pd.to_numeric(raw.iloc[:, 2]), "y": pd.to_numeric(raw.iloc[:, 3]),
                           "kb": pd.to_numeric(raw.iloc[:, 4]), "td": pd.to_numeric(raw.iloc[:, 5])})
        issues.append(Issue("minor", "metadata", "Header file does not state XY or depth units",
                            source=path.name))
    df["xy_unit"], df["depth_unit"] = units["xy_unit"], units["depth_unit"]
    df["source"] = path.name
    return df.reset_index(drop=True), issues


def _dedupe(names: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for n in names:
        if n in seen:
            seen[n] += 1
            out.append(f"{n}_{seen[n]}")
        else:
            seen[n] = 0
            out.append(n)
    return out


def read_deviation(path: str | Path) -> tuple[dict, pd.DataFrame, list[Issue]]:
    """Read a Petrel 'WELL TRACE FROM PETREL' deviation export.

    Returns (metadata, trajectory frame, issues). The trajectory has columns
    MD, X, Y, Z, TVD, DX, DY, AZIM, INCL, DLS (whatever the file provides).
    """
    path = Path(path)
    lines = path.read_text(errors="replace").splitlines()
    meta: dict = {"source": path.name}
    issues: list[Issue] = []
    header_cols: list[str] | None = None
    data = []
    for ln in lines:
        s = ln.strip()
        if s.startswith("#"):
            body = s.lstrip("#").strip()
            if ":" in body and not body.startswith("="):
                k, v = body.split(":", 1)
                k, v = k.strip().upper(), v.strip()
                if k == "WELL NAME":
                    meta["name"] = v
                elif k == "WELL HEAD X-COORDINATE":
                    meta["x"] = float(v.split()[0])
                elif k == "WELL HEAD Y-COORDINATE":
                    meta["y"] = float(v.split()[0])
                elif k == "WELL KB":
                    meta["kb"] = float(v.split()[0])
                    if "(" in v:
                        meta["kb_unit"] = v.split("(")[1].rstrip(")").lower()
                elif k == "WELL TYPE":
                    meta["well_type"] = v
            up = body.upper()
            if "COORDINATE SYSTEM" in up:
                crs = body.split("SYSTEM", 1)[1].strip()
                meta["crs"] = crs
                if crs.strip("<>").lower() in {"undefined", "null", ""}:
                    issues.append(Issue("major", "crs", f"Coordinate reference system is {crs}",
                                        well=meta.get("name"), source=path.name))
            if "AZIMUTH REFERENCE UNDEFINED" in up:
                issues.append(Issue("minor", "crs", "Azimuth reference (true/grid north) undefined",
                                    well=meta.get("name"), source=path.name))
            m = re.search(r"DEPTH \(Z, TVD\) GIVEN IN (\w+)-UNITS", up)
            if m:
                meta["depth_unit"] = m.group(1).lower()
            continue
        if not s:
            continue
        if s.upper().startswith("MD"):
            header_cols = s.split()
            continue
        try:
            data.append([float(v) for v in s.split()])
        except ValueError:
            continue
    if header_cols is None:
        raise ValueError(f"{path.name}: no column header line found")
    traj = pd.DataFrame(data, columns=header_cols[: len(data[0])])
    meta.setdefault("depth_unit", meta.get("kb_unit"))
    if meta.get("depth_unit") is None:
        issues.append(Issue("major", "units", "Deviation file does not state depth units",
                            well=meta.get("name"), source=path.name))
    if not traj["MD"].is_monotonic_increasing:
        issues.append(Issue("critical", "trajectory", "MD is not monotonically increasing",
                            well=meta.get("name"), source=path.name))
    return meta, traj, issues


def read_well_tops(path: str | Path) -> tuple[pd.DataFrame, list[Issue]]:
    """Read a Petrel well-tops export (VERSION 2, BEGIN HEADER ... END HEADER)."""
    path = Path(path)
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    fields, i0 = _parse_header_block(lines)
    rows = [shlex.split(ln) for ln in lines[i0:] if ln.strip()]
    issues: list[Issue] = []
    rows = [r[: len(fields)] + [""] * max(0, len(fields) - len(r)) for r in rows]
    df = pd.DataFrame(rows, columns=_dedupe(fields))
    for c in ["X", "Y", "Z", "MD", "TWT picked", "TWT auto"]:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce").replace({-999: np.nan, -999.25: np.nan})
    units = {m.group(1).lower(): m.group(2).lower() for m in _UNIT_RE.finditer(text)}
    df["depth_unit"] = units.get("depth")
    df["source"] = path.name
    keep = [c for c in ["Well", "Surface", "Type", "MD", "X", "Y", "Z", "TWT picked", "Interpreter",
                        "depth_unit", "source"] if c in df]
    return df[keep].rename(columns={"Well": "well", "Surface": "surface", "Type": "type"}), issues


def read_grid_extent(path: str | Path) -> dict | None:
    """Read the header of a Petrel/IRAP-style ASCII surface (FSASCI) and return extent and stats."""
    path = Path(path)
    lines = path.read_text(errors="replace").splitlines()
    info: dict = {"source": path.name}
    for ln in lines[:20]:
        tok = ln.split()
        if not tok:
            continue
        if tok[0] == "FSLIMI":
            info.update(xmin=float(tok[1]), xmax=float(tok[2]), ymin=float(tok[3]), ymax=float(tok[4]),
                        zmin=float(tok[5]), zmax=float(tok[6]))
        elif tok[0] == "FSNROW":
            info.update(nrow=int(tok[1]), ncol=int(tok[2]))
        elif tok[0] == "FSXINC":
            info.update(xinc=float(tok[1]), yinc=float(tok[2]))
    return info if "xmin" in info else None
