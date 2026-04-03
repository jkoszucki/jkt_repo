"""
Box plot of wGRR similarity to the reference K-locus, per category.

For each locus reads per_locus_dir/{KL}/locus/grr_vs_reference.tsv and
renders a panel with all three categories on the X-axis (match, locus_only, pc_only).
Categories with no data are shown as empty axes.

Category colours:
  match      → green  (#2ca02c)
  locus_only → blue   (#1f77b4)
  pc_only    → yellow (#c9a227)

Outputs:
  plots_dir/grr_vs_reference_all.png     — 18-panel grid, one panel per K-locus
  plots_dir/grr_vs_reference_global.png  — three KDE distributions pooled across all K-loci
"""

from __future__ import annotations

import math
from pathlib import Path

from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CATEGORY_ORDER = ["match", "locus_only", "pc_only"]
CATEGORY_LABELS = {
    "match":      "True positive (TP)",
    "locus_only": "False negative (FN)",
    "pc_only":    "False positive (FP)",
}
CATEGORY_COLORS = {
    "match":      "#2ca02c",
    "locus_only": "#1f77b4",
    "pc_only":    "#c9a227",
}


def _draw_locus_panel(
    ax: plt.Axes,
    locus: str,
    df: pd.DataFrame,
    rng: np.random.Generator,
    tick_fs: int,
    tick_fw: str,
    label_fs: int,
    label_fw: str,
) -> None:
    """Draw a single locus panel (box + jitter) onto *ax*.

    Always draws all three categories on the X-axis; categories without data
    are shown as empty tick positions. Y-axis is always fixed to [0, 1].
    """
    data   = [df.loc[df["category"] == c, "wgrr"].values for c in CATEGORY_ORDER]
    labels = [CATEGORY_LABELS[c] for c in CATEGORY_ORDER]
    colors = [CATEGORY_COLORS[c] for c in CATEGORY_ORDER]

    # Only pass non-empty series to boxplot (positions keep the fixed slots)
    nonempty_idx  = [i for i, d in enumerate(data) if len(d) > 0]
    nonempty_data = [data[i] for i in nonempty_idx]

    if nonempty_data:
        bp = ax.boxplot(
            nonempty_data,
            positions=nonempty_idx,
            patch_artist=True,
            widths=0.5,
            medianprops=dict(color="black", lw=1.5),
            whiskerprops=dict(color="#555555", lw=0.8),
            capprops=dict(color="#555555", lw=0.8),
            showfliers=False,
        )
        for patch, i in zip(bp["boxes"], nonempty_idx):
            patch.set_facecolor(colors[i])
            patch.set_alpha(0.55)

        for i in nonempty_idx:
            x_jitter = rng.uniform(-0.18, 0.18, size=len(data[i])) + i
            ax.scatter(x_jitter, data[i], color=colors[i], alpha=0.55, s=4,
                       zorder=3, linewidths=0)

    ax.set_xlim(-0.5, len(CATEGORY_ORDER) - 0.5)
    ax.set_xticks(range(len(CATEGORY_ORDER)))
    ax.set_xticklabels(labels, fontsize=tick_fs - 1, fontweight=tick_fw, rotation=30, ha="right")
    ax.set_ylim(0, 1)
    ax.set_title(locus, fontsize=label_fs - 1, fontweight=label_fw)
    ax.tick_params(axis="y", labelsize=tick_fs - 1)


def plot_grr_vs_reference_all(
    best_predictors_csv: Path,
    per_locus_dir: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Render all loci in a single grid figure.

    Args:
        best_predictors_csv: cfg.output_dir / "chapter2" / "best_predictors.csv"
        per_locus_dir:       cfg.output_dir / "chapter2" / "per-locus"
        plots_dir:           scripts/figures/chapter2/plots
        style:               cfg.style
    """
    dpi      = getattr(style, "dpi", 300)
    tick_fs  = getattr(style, "tick_fontsize", 8)
    tick_fw  = getattr(style, "tick_fontweight", "bold")
    label_fs = getattr(style, "axis_label_fontsize", 12)
    label_fw = getattr(style, "axis_label_fontweight", "bold")

    best = pd.read_csv(best_predictors_csv)
    rng  = np.random.default_rng(42)

    # Collect loci that have data, then sort by KL number
    loci_data: list[tuple[str, pd.DataFrame]] = []
    for _, row in best.iterrows():
        locus = row["locus"]
        tsv   = per_locus_dir / locus / "locus" / "grr_vs_reference.tsv"
        if not tsv.exists():
            print(f"  [grr_all] {locus}: grr_vs_reference.tsv not found, skipping")
            continue
        loci_data.append((locus, pd.read_csv(tsv, sep="\t")))
    loci_data.sort(key=lambda x: int(x[0].replace("KL", "")))

    n = len(loci_data)
    if n == 0:
        print("  [grr_all] no data found — skipping combined plot.")
        return

    n_cols = 6
    n_rows = math.ceil(n / n_cols)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 2.8, n_rows * 3.2),
        squeeze=False,
    )

    for idx, (locus, df) in enumerate(loci_data):
        ax = axes[idx // n_cols][idx % n_cols]
        _draw_locus_panel(ax, locus, df, rng, tick_fs, tick_fw, label_fs, label_fw)
        if idx % n_cols == 0:
            ax.set_ylabel("wGRR to reference", fontsize=tick_fs, fontweight=tick_fw)

    # Hide unused panels
    for idx in range(n, n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].set_visible(False)

    plt.tight_layout()
    out_path = plots_dir / "grr_vs_reference_all.png"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  [grr_all] → {out_path.name}")


def plot_grr_vs_reference_global(
    best_predictors_csv: Path,
    per_locus_dir: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Pool wGRR values across all K-loci and render three stacked horizontal histograms
    (one per category) sharing the same Y-axis (wGRR 0–1).

    Layout: one column of three subplots (match / locus_only / pc_only).
    X-axis: wGRR to reference K-locus (0–1).
    Y-axis: counts.
    A KDE line is overlaid without shading so histograms remain visible.

    Args:
        best_predictors_csv: cfg.output_dir / "chapter2" / "best_predictors.csv"
        per_locus_dir:       cfg.output_dir / "chapter2" / "per-locus"
        plots_dir:           scripts/figures/chapter2/plots
        style:               cfg.style
    """
    dpi      = getattr(style, "dpi", 300)
    tick_fs  = getattr(style, "tick_fontsize", 8)
    tick_fw  = getattr(style, "tick_fontweight", "bold")
    label_fs = getattr(style, "axis_label_fontsize", 12)
    label_fw = getattr(style, "axis_label_fontweight", "bold")

    best = pd.read_csv(best_predictors_csv)

    # Pool wGRR values across all loci per category
    pooled: dict[str, list[float]] = {c: [] for c in CATEGORY_ORDER}
    for _, row in best.iterrows():
        tsv = per_locus_dir / row["locus"] / "locus" / "grr_vs_reference.tsv"
        if not tsv.exists():
            continue
        df = pd.read_csv(tsv, sep="\t")
        for cat in CATEGORY_ORDER:
            pooled[cat].extend(df.loc[df["category"] == cat, "wgrr"].tolist())

    n_bins  = 50
    y_edges = np.linspace(0, 1, n_bins + 1)
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2

    cat_fs   = round(tick_fs * 0.8)      # category labels — 80% of tick size
    tick_fs  = round(tick_fs  * 0.5)   # tick labels — halved
    label_fs = round(label_fs * 0.5)   # axis labels — halved

    fig, axes = plt.subplots(
        3, 1,
        figsize=(3.75, 2.4),
        sharex=True,
        sharey=True,
        gridspec_kw={"hspace": 0.25},
    )

    for ax, cat in zip(axes, CATEGORY_ORDER):
        vals  = np.array(pooled[cat])
        color = CATEGORY_COLORS[cat]
        label = CATEGORY_LABELS[cat]

        bins   = np.linspace(0, 1, (n_bins // 2) + 1) if cat == "pc_only" else y_edges
        counts, _ = np.histogram(vals, bins=bins)
        centers   = (bins[:-1] + bins[1:]) / 2
        ax.bar(centers, counts, width=(bins[1] - bins[0]) * 0.9,
               color=color, alpha=0.65, edgecolor="white", linewidth=0.4)

        # KDE line (no shading) scaled to counts
        if len(vals) >= 2:
            kde     = gaussian_kde(vals, bw_method="scott")
            y_grid  = np.linspace(0, 1, 400)
            density = kde(y_grid)
            bin_width = y_edges[1] - y_edges[0]
            density_scaled = density * len(vals) * bin_width
            ax.plot(y_grid, density_scaled, color=color, lw=2, alpha=0.75)

        # Median line
        if len(vals) > 0:
            ax.axvline(np.median(vals), color="black", lw=1, linestyle="--", alpha=0.8)

        ax.text(0.02, 0.97, f"{label}\n(n={len(vals)})", transform=ax.transAxes,
                fontsize=cat_fs, fontweight=tick_fw, color="black",
                va="top", ha="left")
        ax.set_xlim(0, 1)
        ax.set_ylabel("")
        ax.tick_params(labelsize=tick_fs)
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_fontweight(tick_fw)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.text(0.02, 0.55, "Count", ha="center", va="center", rotation="vertical",
             fontsize=label_fs, fontweight=label_fw)
    axes[-1].set_xlabel("wGRR to reference K-locus", fontsize=label_fs, fontweight=label_fw)

    fig.subplots_adjust(left=0.10, right=0.97, top=0.93, bottom=0.18)
    out_path = plots_dir / "grr_vs_reference_global.png"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  [grr_global] → {out_path.name}")
