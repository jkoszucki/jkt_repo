from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from pathlib import Path


# Marker per clustering level — mirrors R ggstar shapes
_MARKERS = {
    "PCI80C80": "*",
    "PCI50C50": "D",
    "PCI80C50": "^",
    "PCI50C80": "s",
    "PCI00C80": "P",
    "PCI00C50": "o",
}

_MARKER_LABELS = {
    "PCI80C80": "PCI80C80",
    "PCI50C50": "PCI50C50",
    "PCI80C50": "PCI80C50",
    "PCI50C80": "PCI50C80",
    "PCI00C80": "PCI00C80",
    "PCI00C50": "PCI00C50",
}


def _kl_sort_key(locus: str) -> int:
    try:
        return int(locus.replace("KL", ""))
    except ValueError:
        return 9999


_ALL_K_LOCI = [
    "KL1", "KL2", "KL3", "KL6", "KL7", "KL8", "KL10", "KL11", "KL13", "KL14",
    "KL15", "KL16", "KL21", "KL22", "KL24", "KL25", "KL28", "KL30", "KL35",
    "KL38", "KL39", "KL47", "KL48", "KL55", "KL57", "KL60", "KL61", "KL62",
    "KL64", "KL103", "KL111", "KL122", "KL125", "KL127", "KL142",
]


def plot_sgnh_recall(best_predictors_csv: Path, plots_dir: Path, style=None) -> None:
    df = pd.read_csv(best_predictors_csv)

    # axis label sizes
    axis_label_fs = getattr(style, "axis_label_fontsize", 12)
    axis_label_fw = getattr(style, "axis_label_fontweight", "bold")
    tick_fs       = getattr(style, "tick_fontsize", 8)
    tick_fw       = getattr(style, "tick_fontweight", "bold")
    dpi           = getattr(style, "dpi", 300)
    fill_color    = getattr(style, "sgnh_domain_color", "#c9a227")

    x_labels = sorted(_ALL_K_LOCI, key=_kl_sort_key)
    x_positions = {locus: i for i, locus in enumerate(x_labels)}

    # 4 discrete size classes for PC_abundance
    _SIZE_BINS   = [(1, 10), (11, 25), (26, float("inf"))]
    _SIZE_VALUES = [50, 150, 400]
    _SIZE_LABELS = ["1–10", "11–25", ">25"]

    def _abundance_to_size(ab):
        for (lo, hi), s in zip(_SIZE_BINS, _SIZE_VALUES):
            if lo <= ab <= hi:
                return s
        return _SIZE_VALUES[-1]

    fig, ax = plt.subplots(figsize=(10, 4))

    for _, row in df.iterrows():
        marker = _MARKERS.get(row["version"], "o")
        ax.scatter(
            x_positions[row["locus"]], row["recall"],
            marker=marker,
            s=_abundance_to_size(row["PC_abundance"]),
            color=fill_color,
            edgecolors="black",
            linewidths=0.8,
            zorder=3,
            alpha=0.85,
        )

    ax.set_xlim(-0.5, len(x_labels) - 0.5)
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels)
    ax.set_ylim(0, 1)
    ax.set_xlabel("K-locus", fontsize=axis_label_fs, fontweight=axis_label_fw)
    ax.set_ylabel("Recall", fontsize=axis_label_fs, fontweight=axis_label_fw)
    ax.tick_params(axis="x", labelsize=tick_fs, labelrotation=45)
    ax.tick_params(axis="y", labelsize=tick_fs)
    loci_with_hits = set(df["locus"])
    for label in ax.get_xticklabels():
        if label.get_text() in loci_with_hits:
            label.set_fontweight(tick_fw)
        else:
            label.set_fontweight("light")
            label.set_color("gray")
    for label in ax.get_yticklabels():
        label.set_fontweight(tick_fw)

    ax.grid(axis="both", linestyle="--", linewidth=0.5, alpha=0.4, color="gray")
    ax.set_axisbelow(True)

    # legend: clustering level → marker shape
    versions_present = df["version"].unique()
    shape_handles = [
        mpatches.Patch(visible=False)  # spacer
    ]
    shape_handles = []
    for v in sorted(versions_present):
        m = _MARKERS.get(v, "o")
        h = ax.scatter([], [], marker=m, s=80, color=fill_color,
                       edgecolors="black", linewidths=0.8, label=v)
        shape_handles.append(h)

    shape_legend = ax.legend(
        handles=shape_handles,
        title="Clustering level",
        title_fontsize=tick_fs,
        fontsize=tick_fs,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        borderaxespad=0,
        frameon=True,
    )
    ax.add_artist(shape_legend)

    # Size legend: 4 discrete classes
    size_handles = [
        ax.scatter([], [], marker="o", s=s,
                   color=fill_color, edgecolors="black", linewidths=0.8,
                   label=label)
        for s, label in zip(_SIZE_VALUES, _SIZE_LABELS)
    ]
    ax.legend(
        handles=size_handles,
        title="Number of sequences",
        title_fontsize=tick_fs,
        fontsize=tick_fs,
        loc="lower left",
        bbox_to_anchor=(1.01, 0),
        borderaxespad=0,
        frameon=True,
        labelspacing=1.5,
        handletextpad=1.0,
    )

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = plots_dir / f"figure_sgnh_recall.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"Saved → {out}")
    plt.close(fig)
