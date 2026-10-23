"""Static QC figures (matplotlib). Colours follow a validated colour-blind-safe categorical order."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
CRITICAL = "#d03b3b"
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#a9a8a3", "#e4e3df", "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.titlecolor": INK, "font.size": 9, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
})


def log_panel(raw: pd.DataFrame, clean: pd.DataFrame, well: str, unit: str, path: Path,
              curves=("GR", "DT", "RHOB", "NPHI"), pred: dict[str, pd.Series] | None = None,
              interval: tuple[float, float] | None = None) -> Path:
    curves = [c for c in curves if c in raw.columns or (pred and c in pred)]
    fig, axes = plt.subplots(1, len(curves), figsize=(2.3 * len(curves), 7.5), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, c in zip(axes, curves, strict=True):
        if c in raw:
            r = raw[c]
            if interval:
                r = r.loc[interval[0]:interval[1]]
            ax.plot(r.to_numpy(), r.index, color=MUTED, lw=0.8, label="raw")
        if c in clean:
            s = clean[c]
            if interval:
                s = s.loc[interval[0]:interval[1]]
            ax.plot(s.to_numpy(), s.index, color=BLUE, lw=1.0, label="conditioned")
        if pred and c in pred and pred[c] is not None:
            p = pred[c]
            if c in clean:
                p = p.where(clean[c].reindex(p.index).isna())  # show prediction only where data is missing
            if interval:
                p = p.loc[interval[0]:interval[1]]
            if p.notna().any():
                ax.plot(p.to_numpy(), p.index, color=ORANGE, lw=1.0, label="predicted (gap fill)")
        ref = clean[c] if c in clean else (pred or {}).get(c)
        if ref is not None and ref.notna().any():
            r = ref.loc[interval[0]:interval[1]] if interval else ref
            lo, hi = np.nanpercentile(r, [0.5, 99.5])
            pad = 0.15 * (hi - lo)
            ax.set_xlim(lo - pad, hi + pad)
        ax.set_title(c)
        ax.xaxis.set_label_position("top")
        if c == "DT":
            ax.invert_xaxis()
    axes[0].invert_yaxis()
    axes[0].set_ylabel(f"MD ({unit})")
    handles, labels = [], []
    for ax in axes:
        for h, lab in zip(*ax.get_legend_handles_labels(), strict=True):
            if lab not in labels:
                handles.append(h)
                labels.append(lab)
    fig.legend(handles, labels, loc="lower center", ncol=len(labels), bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(f"{well}: raw vs conditioned logs", fontweight="bold", color=INK, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def rpt_crossplot(samples: pd.DataFrame, template: dict[str, pd.DataFrame], path: Path,
                  title: str = "Rock physics template QC", flag_col: str | None = "RPT_FLAG") -> Path:
    fig, ax = plt.subplots(figsize=(7, 5.2))
    ok = ~samples[flag_col].astype(bool) if flag_col in samples else np.ones(len(samples), bool)
    ax.scatter(samples.loc[ok, "AI"], samples.loc[ok, "VpVs"], s=4, color=MUTED, alpha=0.35, lw=0,
               label="samples on a template trend")
    if flag_col in samples and (~ok).any():
        ax.scatter(samples.loc[~ok, "AI"], samples.loc[~ok, "VpVs"], s=10, color=CRITICAL, lw=0.8,
                   marker="x", label=f"off-template (flagged): {int((~ok).sum())}")
    for (name, df), col in zip(template.items(), (BLUE, ORANGE, AQUA), strict=False):
        ax.plot(df.AI, df.VpVs, color=col, lw=2, label=name)
        ax.annotate(name, (df.AI.iloc[-1], df.VpVs.iloc[-1]), color=INK2, fontsize=8,
                    xytext=(4, 0), textcoords="offset points", va="center")
    ax.set_xlabel("Acoustic impedance (km/s x g/cm3)")
    ax.set_ylabel("Vp/Vs")
    ax.set_title(title, loc="left")
    ax.legend(loc="upper right", markerscale=2)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def blind_well(measured: pd.Series, predicted: pd.Series, well: str, r2: float, unit: str, path: Path,
               curve: str = "DT") -> Path:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.5, 6.5), gridspec_kw={"width_ratios": [1, 1.3]})
    ax1.plot(measured.to_numpy(), measured.index, color=MUTED, lw=0.9, label="measured")
    ax1.plot(predicted.to_numpy(), predicted.index, color=BLUE, lw=0.9, label="predicted (blind)")
    ax1.invert_yaxis()
    ax1.invert_xaxis()
    ax1.set_ylabel(f"MD ({unit})")
    ax1.set_title(f"{curve} track", loc="left")
    ax1.legend(loc="lower left")
    j = pd.concat([measured, predicted], axis=1).dropna()
    ax2.scatter(j.iloc[:, 0], j.iloc[:, 1], s=3, color=BLUE, alpha=0.3, lw=0)
    lo, hi = np.percentile(j.to_numpy(), [0.5, 99.5])
    ax2.plot([lo, hi], [lo, hi], color=INK2, lw=1, ls="--")
    ax2.annotate("1:1", (hi, hi), color=INK2, fontsize=8, xytext=(-18, -10), textcoords="offset points")
    ax2.set_xlabel(f"measured {curve} (us/ft)")
    ax2.set_ylabel(f"predicted {curve} (us/ft)")
    ax2.set_title(f"R² = {r2:.2f}", loc="left")
    fig.suptitle(f"Blind-well test: {well} held out of training", fontweight="bold", color=INK, x=0.02,
                 ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def checkshot_qc(df: pd.DataFrame, well: str, unit: str, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(5, 6.5))
    flag = df["FLAG"]
    # VINT[i] spans levels i-1 -> i, so a mis-picked level i affects VINT[i] and VINT[i + 1]
    mis = (flag == "mis-picked level") | (flag.shift(1) == "mis-picked level")
    ano = flag == "anomalous interval"
    ok = ~(mis | ano)
    mid = (df["DEPTH"] + df["DEPTH"].shift(1)) / 2  # plot each interval at its mid-depth
    ax.scatter(df.loc[ok, "VINT"], mid[ok], s=14, color=MUTED, label="interval velocity")
    if mis.any():
        ax.scatter(df.loc[mis, "VINT"], mid[mis], s=34, color=CRITICAL, marker="x", lw=1.2,
                   label="mis-picked level (both intervals)")
    if ano.any():
        ax.scatter(df.loc[ano, "VINT"], mid[ano], s=30, facecolors="none", edgecolors=ORANGE, lw=1.2,
                   label="anomalous interval (verify)")
    ax.plot(df["VINT_SMOOTH"], df["DEPTH"], color=BLUE, lw=2, label="smoothed (quantisation-aware)")
    ax.invert_yaxis()
    ax.set_xlabel(f"Interval velocity ({unit}/s)")
    ax.set_ylabel(f"Depth ({unit})")
    ax.set_title(f"{well}: checkshot QC", loc="left")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def basemap(master: pd.DataFrame, surveys: list[dict], path: Path, duplicates: set[str] | None = None,
            title: str = "Basemap") -> Path:
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for s, col in zip(surveys, (BLUE, ORANGE, AQUA), strict=False):
        g = s["grid"]
        from .io.segy import grid_to_xy
        il = [s["il"][0], s["il"][1], s["il"][1], s["il"][0], s["il"][0]]
        xl = [s["xl"][0], s["xl"][0], s["xl"][1], s["xl"][1], s["xl"][0]]
        x, y = grid_to_xy(g, il, xl)
        ax.plot(x, y, color=col, lw=2, label=f"3D survey: {s['name']}")
    placed: list[tuple[float, float]] = []
    span = max(np.ptp(master.x.to_numpy(float)), np.ptp(master.y.to_numpy(float)), 1.0)
    for _, r in master.sort_values("x").iterrows():
        isdup = r["n_sources"] > 1
        ax.scatter(r.x, r.y, s=40 if isdup else 24, color=CRITICAL if isdup else INK, zorder=3,
                   marker="D" if isdup else "o")
        near = sum(np.hypot(r.x - px, r.y - py) < 0.15 * span for px, py in placed)
        dy = 6 if near % 2 == 0 else -12
        ax.annotate(r["standard_name"] + (" (duplicate)" if isdup else ""), (r.x, r.y), xytext=(6, dy),
                    textcoords="offset points", fontsize=7.5, color=INK2)
        placed.append((r.x, r.y))
    ax.set_aspect("equal", adjustable="datalim")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_title(title, loc="left")
    if surveys:
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=len(surveys))
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path
