"""
Chapter 2, Figure 2, Panel B — precision vs recall scatter, 6-panel grid.

One subplot per MMseqs2 clustering level.
Data: pyseer_hits_sslbh_pvalcor005.tsv (lasso, pvalue_corr ≤ 0.05).

Colour logic:
  acetyltransferase (acetyl_hhsearch or enzyme_tmscore_05) → sslbh_color, alpha=0.75
  rest                                                      → gray_color,  alpha=0.35
  acetyltransferase with enzyme_tmscore_05=True             → thick colored outline
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


_CLUSTERING_LEVELS = [
    ("PCI00C50", "0% id / 50% cov"),
    ("PCI50C50", "50% id / 50% cov"),
    ("PCI80C50", "80% id / 50% cov"),
    ("PCI00C80", "0% id / 80% cov"),
    ("PCI50C80", "50% id / 80% cov"),
    ("PCI80C80", "80% id / 80% cov"),
]


def plot_figure2_panelB(
    pyseer_tsv: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Produce Figure 2 Panel B: 6-panel precision vs recall scatter.

    Args:
        pyseer_tsv: pyseer_hits_sslbh_pvalcor005.tsv
        plots_dir:  output directory for plots
        style:      cfg.style (optional)
    """
    dpi      = getattr(style, "dpi", 300)
    tick_fs  = getattr(style, "tick_fontsize", 8)
    tick_fw  = getattr(style, "tick_fontweight", "bold")
    label_fs = getattr(style, "axis_label_fontsize", 12)
    label_fw = getattr(style, "axis_label_fontweight", "bold")

    sslbh_color = getattr(style, "sslbh_color", "#9467bd")
    gray_color  = getattr(style, "gray_color",  "#bfbfbf")

    df = pd.read_csv(pyseer_tsv, sep="\t")

    is_acetyl   = df["acetyl_hhsearch"] | df["enzyme_tmscore_05"]
    is_struct   = df["enzyme_tmscore_05"]

    df["_color"]  = df.apply(lambda r: sslbh_color if is_acetyl[r.name] else gray_color, axis=1)
    df["_alpha"]  = df.apply(lambda r: 0.75 if is_acetyl[r.name] else 0.35, axis=1)
    df["_struct"] = is_struct

    fig, axes = plt.subplots(2, 3, figsize=(11, 7.5), sharex=True, sharey=True)
    axes_flat = axes.flatten()

    for ax_idx, (cl, label) in enumerate(_CLUSTERING_LEVELS):
        ax  = axes_flat[ax_idx]
        sub = df[df["clustering_level"] == cl]

        bg   = sub[~is_acetyl[sub.index]]
        ac   = sub[is_acetyl[sub.index] & ~is_struct[sub.index]]
        strc = sub[is_struct[sub.index]]

        ax.scatter(bg["precision"],   bg["recall"],
                   c=gray_color,  s=60,  alpha=0.35,
                   edgecolors="none", zorder=1)
        ax.scatter(ac["precision"],   ac["recall"],
                   c=sslbh_color, s=100, alpha=0.75,
                   edgecolors="none", zorder=2)
        ax.scatter(strc["precision"], strc["recall"],
                   c=sslbh_color, s=120, alpha=0.75,
                   edgecolors="black", linewidths=2.0,
                   zorder=3)

        n_pc    = sub["PC"].nunique()
        n_acetyl = (is_acetyl[sub.index]).sum()

        ax.set_title(label, fontsize=tick_fs, fontweight=tick_fw)
        ax.text(0.97, 0.03,
                f"n={n_pc} PCs\nacetyl={n_acetyl}",
                ha="right", va="bottom", transform=ax.transAxes,
                fontsize=tick_fs - 1, color="#555555")

        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.03)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=tick_fs - 1)
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_fontweight(tick_fw)

    fig.supxlabel("Precision", fontsize=round(label_fs * 0.85), fontweight=label_fw)
    fig.supylabel("Recall",    fontsize=round(label_fs * 0.85), fontweight=label_fw)

    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    def _save_legend(handles, stem):
        fig_leg, ax_leg = plt.subplots(figsize=(3, len(handles) * 0.35 + 0.2))
        ax_leg.axis("off")
        ax_leg.legend(handles=handles, loc="center", ncol=1,
                      fontsize=tick_fs, frameon=False)
        for ext in ("png", "pdf"):
            out = plots_dir / f"{stem}.{ext}"
            fig_leg.savefig(out, dpi=dpi, bbox_inches="tight")
            print(f"  [{stem}] → {out.name}")
        plt.close(fig_leg)
    legend_elements = [
        mpatches.Patch(facecolor=sslbh_color, edgecolor="none",
                       label="Acetyltransferase (HHsearch)"),
        mlines.Line2D([], [], marker="o", color="none",
                      markerfacecolor=sslbh_color,
                      markeredgecolor="black", markeredgewidth=2.0,
                      markersize=8, label="+ structural similarity (TM-score ≥ 0.5)"),
        mpatches.Patch(facecolor=gray_color, edgecolor="none",
                       alpha=0.35, label="Other"),
    ]
    plots_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        out = plots_dir / f"figure2-panelB.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure2-panelB] → {out.name}")
    plt.close(fig)

    _save_legend(legend_elements, "figure2-panelB-legend")
