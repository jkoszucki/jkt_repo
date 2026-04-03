"""
Chapter 2, Figure 1, Panel A — 6-panel F1 vs MCC scatter.

One subplot per MMseqs2 clustering level, coloured by ECOD topology.
Input: gwas_hits.tsv (pre-filtered and pre-classified by processing/gwas-proc).

Layout: 3 columns × 2 rows
  Row 1:  PCI00C50  |  PCI50C50  |  PCI80C50
  Row 2:  PCI00C80  |  PCI50C80  |  PCI80C80

Colour mapping (from ecod_type + reported_topology_PC):
  sgnh-ecod                              → sgnh_domain_color  (#c9a227)
  ssrbh-ecod  (Pectin lyase-like)        → ssrbh_color        (#1f77b4)
  other-ecod  with SSLBH topology        → sslbh_color        (#9467bd)
  other-ecod  any other topology         → gray_color         (#bfbfbf)
  no-ecod                                → no_ecod_color      (#ffffff, gray border)
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd


def _save_legend(handles, plots_dir, stem, fontsize, dpi):
    """Save legend elements as a standalone figure."""
    fig_leg, ax_leg = plt.subplots(figsize=(3, len(handles) * 0.35 + 0.2))
    ax_leg.axis("off")
    ax_leg.legend(handles=handles, loc="center", ncol=1,
                  fontsize=fontsize, frameon=False)
    for ext in ("png", "pdf"):
        out = plots_dir / f"{stem}.{ext}"
        fig_leg.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [{stem}] → {out.name}")
    plt.close(fig_leg)


_CLUSTERING_LEVELS = [
    ("PCI00C50", "0% id / 50% cov"),
    ("PCI50C50", "50% id / 50% cov"),
    ("PCI80C50", "80% id / 50% cov"),
    ("PCI00C80", "0% id / 80% cov"),
    ("PCI50C80", "50% id / 80% cov"),
    ("PCI80C80", "80% id / 80% cov"),
]


def _ecod_color(ecod_type: str, reported_topology: str, style) -> str:
    sgnh_color    = getattr(style, "sgnh_domain_color", "#c9a227")
    ssrbh_color   = getattr(style, "ssrbh_color",       "#1f77b4")
    sslbh_color   = getattr(style, "sslbh_color",       "#9467bd")
    gray_color    = getattr(style, "gray_color",         "#bfbfbf")
    no_ecod_color = getattr(style, "no_ecod_color",      "#ffffff")

    if ecod_type == "sgnh-ecod":
        return sgnh_color
    if ecod_type == "ssrbh-ecod":
        return ssrbh_color
    if ecod_type == "no-ecod":
        return no_ecod_color
    # other-ecod: check for SSLBH
    t = str(reported_topology).lower()
    if "left-handed" in t or "sslbh" in t:
        return sslbh_color
    return gray_color


def plot_figure1_panelA(
    gwas_hits_tsv: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Produce Figure 1 Panel A: 6-panel F1 vs MCC scatter coloured by ECOD topology.

    Args:
        gwas_hits_tsv: path to gwas_hits.tsv (output of processing/gwas-proc)
        plots_dir:     output directory for plots
        style:         cfg.style (optional)
    """
    dpi      = getattr(style, "dpi", 300)
    tick_fs  = getattr(style, "tick_fontsize", 8)
    tick_fw  = getattr(style, "tick_fontweight", "bold")
    label_fs = getattr(style, "axis_label_fontsize", 12)
    label_fw = getattr(style, "axis_label_fontweight", "bold")

    sgnh_color    = getattr(style, "sgnh_domain_color", "#c9a227")
    ssrbh_color   = getattr(style, "ssrbh_color",       "#1f77b4")
    sslbh_color   = getattr(style, "sslbh_color",       "#9467bd")
    gray_color    = getattr(style, "gray_color",         "#bfbfbf")
    no_ecod_color = getattr(style, "no_ecod_color",      "#ffffff")

    draw_order = {
        gray_color:    1,
        ssrbh_color:   2,
        sslbh_color:   3,
        sgnh_color:    4,
        no_ecod_color: 5,
    }

    hits = pd.read_csv(gwas_hits_tsv, sep="\t")
    hits["color"] = hits.apply(
        lambda r: _ecod_color(r["ecod_type"], r["reported_topology_PC"], style),
        axis=1,
    )

    fig, axes = plt.subplots(2, 3, figsize=(11, 7.5), sharex=True, sharey=True)
    axes_flat = axes.flatten()

    for ax_idx, (cl, label) in enumerate(_CLUSTERING_LEVELS):
        ax = axes_flat[ax_idx]
        df = hits[hits["clustering_level"] == cl]

        ax.plot([0.47, 1.03], [0.47, 1.03], color="#aaaaaa", lw=0.6,
                linestyle="--", zorder=0)

        for color, subset in sorted(
            df.groupby("color", sort=False),
            key=lambda kv: draw_order.get(kv[0], 2),
        ):
            is_no_ecod = color == no_ecod_color
            ax.scatter(
                subset["F1_score"], subset["MCC"],
                c=color, s=175, alpha=0.75,
                zorder=draw_order.get(color, 2),
                edgecolors="#888888" if is_no_ecod else "white",
                linewidths=0.5 if is_no_ecod else 0.3,
            )

        n_pc      = df["PC"].nunique()
        n_loci    = df["locus"].nunique()
        median_f1 = df["F1_score"].median()

        ax.axvline(median_f1, color="#555555", lw=0.8, linestyle=":", zorder=1)

        ax.set_title(label, fontsize=tick_fs, fontweight=tick_fw)
        ax.text(0.97, 0.03, f"n={n_pc} PCs\n{n_loci} loci",
                ha="right", va="bottom", transform=ax.transAxes,
                fontsize=tick_fs - 1, color="#555555")
        ax.text(median_f1, 0.9, f"med={median_f1:.2f}",
                ha="right", va="center", rotation=90,
                fontsize=tick_fs - 2, color="#555555")

        ax.set_xlim(0.47, 1.03)
        ax.set_ylim(0.47, 1.03)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=tick_fs - 1)
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_fontweight(tick_fw)

    fig.supxlabel("F1 score", fontsize=round(label_fs * 0.85), fontweight=label_fw)
    fig.supylabel("MCC", fontsize=round(label_fs * 0.85), fontweight=label_fw)

    legend_elements = [
        mpatches.Patch(facecolor=sgnh_color,    edgecolor="white", label="SGNH hydrolase"),
        mpatches.Patch(facecolor=ssrbh_color,   edgecolor="white", label="SSRBH (depolymerase)"),
        mpatches.Patch(facecolor=sslbh_color,   edgecolor="white", label="SSLBH (acetyltransferase)"),
        mpatches.Patch(facecolor=gray_color,    edgecolor="white", label="Other ECOD"),
        mpatches.Patch(facecolor=no_ecod_color, edgecolor="#888888",
                       linewidth=0.6, label="No ECOD hit"),
    ]
    plots_dir.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        out = plots_dir / f"figure1-panelA.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure1-panelA] → {out.name}")
    plt.close(fig)

    # Save legend separately
    _save_legend(legend_elements, plots_dir, "figure1-panelA-legend", tick_fs, dpi)
