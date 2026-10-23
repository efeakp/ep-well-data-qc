"""Core data structures: wells, issues and processing lineage."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone

import pandas as pd


@dataclass
class Issue:
    """A single data-quality finding."""

    severity: str  # "critical" | "major" | "minor" | "info"
    category: str  # e.g. "naming", "units", "spike", "duplicate", "segy"
    message: str
    well: str | None = None
    curve: str | None = None
    source: str | None = None

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class LineageStep:
    step: str
    action: str
    well: str | None = None
    curve: str | None = None
    n_affected: int | None = None
    params: str | None = None
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds"))


class Lineage:
    """Append-only record of every processing step, for audit and documentation."""

    def __init__(self) -> None:
        self.steps: list[LineageStep] = []

    def add(self, step: str, action: str, well: str | None = None, curve: str | None = None,
            n_affected: int | None = None, **params) -> None:
        p = ", ".join(f"{k}={v}" for k, v in params.items()) or None
        self.steps.append(LineageStep(step, action, well, curve, n_affected, p))

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame([asdict(s) for s in self.steps])


@dataclass
class Well:
    """A well assembled from one or more source files.

    ``logs`` is indexed by measured depth (MD) in ``depth_unit``; column names are canonical
    mnemonics once :func:`wellqc.standardise.standardise_curves` has run.
    """

    name: str
    uwi: str | None = None
    x: float | None = None
    y: float | None = None
    kb: float | None = None
    td: float | None = None
    depth_unit: str = "m"
    project: str | None = None
    logs: pd.DataFrame | None = None
    units: dict[str, str] = field(default_factory=dict)
    trajectory: pd.DataFrame | None = None  # columns MD, X, Y, Z, TVD, INCL, AZIM
    tops: pd.DataFrame | None = None
    checkshot: pd.DataFrame | None = None  # columns DEPTH, TWT
    aliases: set[str] = field(default_factory=set)
    sources: list[str] = field(default_factory=list)
    master_uwi: str | None = None

    def md_to_tvdss(self, md):
        """Interpolate TVDSS (positive down, below datum) from the trajectory, or assume vertical."""
        import numpy as np

        md = np.asarray(md, dtype=float)
        if self.trajectory is not None and len(self.trajectory) > 1:
            return -np.interp(md, self.trajectory["MD"], self.trajectory["Z"])
        return md - (self.kb or 0.0)
