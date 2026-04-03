"""
Chapter 2, Figure 1, Panel A
─────────────────────────────
6-panel scatter (F1 vs MCC) — one subplot per MMseqs2 clustering level.

Layout: 3 columns × 2 rows
  Row 1:  PCI00C50  |  PCI50C50  |  PCI80C50
  Row 2:  PCI00C80  |  PCI50C80  |  PCI80C80

Each subplot shows all pyseer hits for that clustering level with
F1 ≥ 0.5 and MCC ≥ 0.5, coloured by ECOD topology (reported_topology_PC).

Topology colour mapping (from cfg.style, with fallbacks):
  SGNH hydrolase                          → sgnh_domain_color  (#c9a227)
  Pectin lyase-like (= SSRBH fold)        → ssrbh_color        (#1f77b4)
  single-stranded right-handed beta helix → ssrbh_color        (#1f77b4)
  single-stranded left-handed beta helix  → sslbh_color        (#9467bd)
  No ECOD hit (prob < 0.7 or no entry)   → no_ecod_color      (#ffffff, gray border)
  Any other topology                      → gray_color          (#bfbfbf)

Outputs:
  plots_dir/figure1-panelA.png
  plots_dir/figure1-panelA.pdf
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd


# ---------------------------------------------------------------------------
# Clustering level definitions
# ---------------------------------------------------------------------------

_CLUSTERING_LEVELS = [
    ("PCI00C50", "0% id / 50% cov"),
    ("PCI50C50", "50% id / 50% cov"),
    ("PCI80C50", "80% id / 50% cov"),
    ("PCI00C80", "0% id / 80% cov"),
    ("PCI50C80", "50% id / 80% cov"),
    ("PCI80C80", "80% id / 80% cov"),
]

# Minimum thresholds for the HHsearch annotation to count as a topology hit
_ECOD_PROB_MIN = 0.7
_ECOD_TCOV_MIN = 0.1


# ---------------------------------------------------------------------------
# Annotation helpers
# ---------------------------------------------------------------------------

def _load_annotations(hhsearch_dir: Path, version: str) -> pd.DataFrame:
    """Load clusters_functions.tsv for a given clustering level."""
    tsv = hhsearch_dir / version / "clusters_functions.tsv"
    if not tsv.exists():
        return pd.DataFrame(columns=["PC", "db", "reported_topology_PC", "prob", "tcov"])
    return pd.read_csv(tsv, sep="\t")


def _build_pc_topology_map(ann: pd.DataFrame) -> dict[str, str | None]:
    """
    Return {PC: reported_topology_PC} for PCs with an ECOD hit passing
    the probability and target-coverage thresholds.  PCs with no qualifying
    ECOD hit map to None (= no-ecod).
    """
    ecod = ann[
        (ann["db"] == "ECOD") &
        (ann["prob"] >= _ECOD_PROB_MIN) &
        (ann["tcov"] >= _ECOD_TCOV_MIN)
    ].copy()

    # For PCs with multiple qualifying hits, keep the one with highest prob.
    ecod = ecod.sort_values("prob", ascending=False).drop_duplicates(subset="PC")
    return dict(zip(ecod["PC"], ecod["reported_topology_PC"]))


def _topology_to_color(topology: str | None, style) -> str:
    """
    Map a reported_topology_PC string (or None) to a fill colour.

    "Pectin-lyase like" is the ECOD annotation for the single-stranded
    right-handed beta-helix (SSRBH) fold — the tail fiber / depolymerase
    fold.  It is coloured with ssrbh_color.
    """
    if topology is None:
        return getattr(style, "no_ecod_color", "#ffffff")

    t = str(topology).lower()
    if "sgnh" in t:
        return getattr(style, "sgnh_domain_color", "#c9a227")
    if "pectin" in t:                          # Pectin-lyase like = SSRBH fold
        return getattr(style, "ssrbh_color", "#1f77b4")
    if "left-handed" in t or "sslbh" in t:
        return getattr(style, "sslbh_color", "#9467bd")
    return getattr(style, "gray_color", "#bfbfbf")


# ---------------------------------------------------------------------------
# Panel A — 6-panel scatter
# ---------------------------------------------------------------------------

def plot_gwas_scatter(
    pyseer_hits_tsv: Path,
    hhsearch_dir: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Produce Figure 1 Panel A: 6-panel F1 vs MCC scatter.

    Args:
        pyseer_hits_tsv:  path to pyseer_hits_all.tsv
        hhsearch_dir:     directory containing per-version HHsearch annotations
        plots_dir:        output directory for plots
        style:            cfg.style (optional)
    """
    dpi      = getattr(style, "dpi", 300)
    tick_fs  = getattr(style, "tick_fontsize", 8)
    tick_fw  = getattr(style, "tick_fontweight", "bold")
    label_fs = getattr(style, "axis_label_fontsize", 12)
    label_fw = getattr(style, "axis_label_fontweight", "bold")

    no_ecod_color   = getattr(style, "no_ecod_color",    "#ffffff")
    sgnh_color      = getattr(style, "sgnh_domain_color","#c9a227")
    pectin_color    = getattr(style, "ssrbh_color",       "#1f77b4")
    ssrbh_color     = getattr(style, "ssrbh_color",      "#1f77b4")
    sslbh_color     = getattr(style, "sslbh_color",      "#9467bd")
    gray_color      = getattr(style, "gray_color",        "#bfbfbf")

    hits_all = pd.read_csv(pyseer_hits_tsv, sep="\t")

    fig, axes = plt.subplots(2, 3, figsize=(11, 7.5), sharex=True, sharey=True)
    axes_flat = axes.flatten()

    for ax_idx, (version, label) in enumerate(_CLUSTERING_LEVELS):
        ax = axes_flat[ax_idx]

        df = hits_all[
            (hits_all["version"] == version) &
            (hits_all["F1_score"] >= 0.5) &
            (hits_all["MCC"] >= 0.5)
        ].copy()

        ann       = _load_annotations(hhsearch_dir, version)
        topo_map  = _build_pc_topology_map(ann)

        df["topology"] = df["PC"].map(topo_map)  # None where no ECOD hit
        df["color"]    = df["topology"].apply(lambda t: _topology_to_color(t, style))

        # Draw in order: gray → other → SGNH/SSRBH/SSLBH on top; no-ecod last (white)
        draw_order = {gray_color: 1, pectin_color: 2, ssrbh_color: 3,
                      sslbh_color: 3, sgnh_color: 4, no_ecod_color: 5}

        ax.plot([0.47, 1.03], [0.47, 1.03], color="#aaaaaa", lw=0.6,
                linestyle="--", zorder=0)

        for color, subset in sorted(df.groupby("color", sort=False),
                                    key=lambda kv: draw_order.get(kv[0], 2)):
            is_no_ecod = color == no_ecod_color
            ax.scatter(
                subset["F1_score"], subset["MCC"],
                c=color, s=35, alpha=0.75,
                zorder=draw_order.get(color, 2),
                edgecolors="#888888" if is_no_ecod else "white",
                linewidths=0.5 if is_no_ecod else 0.3,
            )

        n_pc   = df["PC"].nunique()
        n_loci = df["locus"].nunique()
        ax.set_title(label, fontsize=tick_fs, fontweight=tick_fw)
        ax.text(0.97, 0.03, f"n={n_pc} PCs\n{n_loci} loci",
                ha="right", va="bottom", transform=ax.transAxes,
                fontsize=tick_fs - 1, color="#555555")

        ax.set_xlim(0.47, 1.03)
        ax.set_ylim(0.47, 1.03)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=tick_fs - 1)
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_fontweight(tick_fw)

    # Shared axis labels on outer subplots
    for ax in axes[1]:   # bottom row
        ax.set_xlabel("F1 score", fontsize=round(label_fs * 0.85), fontweight=label_fw)
    for ax in axes[:, 0]:  # left column
        ax.set_ylabel("MCC", fontsize=round(label_fs * 0.85), fontweight=label_fw)

    # Legend (one shared legend on the figure)
    legend_elements = [
        mpatches.Patch(facecolor=sgnh_color,    edgecolor="white", label="SGNH hydrolase"),
        mpatches.Patch(facecolor=ssrbh_color,   edgecolor="white",
                       label="SSRBH (depolymerase)"),
        mpatches.Patch(facecolor=sslbh_color,   edgecolor="white", label="SSLBH (acetyltransferase)"),
        mpatches.Patch(facecolor=gray_color,    edgecolor="white", label="Other ECOD"),
        mpatches.Patch(facecolor=no_ecod_color, edgecolor="#888888",
                       linewidth=0.6, label="No ECOD hit"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=3,
        fontsize=tick_fs,
        frameon=False,
        bbox_to_anchor=(0.5, -0.04),
    )

    plt.tight_layout(rect=[0, 0.06, 1, 1])
    for ext in ("png", "pdf"):
        out = plots_dir / f"figure1-panelA.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure1-panelA] → {out.name}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# PHROG functions bar plot (kept for reference / exploratory use)
# ---------------------------------------------------------------------------

def plot_gwas_phrog_barplot(
    pyseer_hits_tsv: Path,
    hhsearch_dir: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Horizontal bar plot of PHROG function counts for PCI50C50 filtered set
    (F1 ≥ 0.5, MCC ≥ 0.5), sorted by frequency.

    Output: plots_dir/gwas_phrog_functions.png
    """
    _VERSION = "PCI50C50"

    dpi      = getattr(style, "dpi", 300)
    tick_fs  = getattr(style, "tick_fontsize", 8)
    tick_fw  = getattr(style, "tick_fontweight", "bold")
    label_fs = getattr(style, "axis_label_fontsize", 12)
    label_fw = getattr(style, "axis_label_fontweight", "bold")
    gray_color = getattr(style, "gray_color", "#bfbfbf")

    hits = pd.read_csv(pyseer_hits_tsv, sep="\t")
    pcs = hits[
        (hits["version"] == _VERSION) &
        (hits["F1_score"] >= 0.5) &
        (hits["MCC"] >= 0.5)
    ]["PC"].unique()

    ann    = _load_annotations(hhsearch_dir, _VERSION)
    phrogs = ann[ann["db"] == "PHROGS"].copy()

    func_counts: dict[str, int] = {}
    for pc in pcs:
        sub = phrogs[phrogs["PC"] == pc]
        func = str(sub.loc[sub["bits"].idxmax(), "function"]) if not sub.empty else "no PHROG hit"
        func_counts[func] = func_counts.get(func, 0) + 1

    funcs  = sorted(func_counts, key=func_counts.get, reverse=True)
    counts = [func_counts[f] for f in funcs]

    fig, ax = plt.subplots(figsize=(4.0, max(2.5, len(funcs) * 0.6)))
    ax.barh(funcs, counts, color=gray_color, edgecolor="white", linewidth=0.4)

    ax.set_xlabel("Number of PCs", fontsize=round(label_fs * 0.8), fontweight=label_fw)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.tick_params(axis="y", labelsize=tick_fs)
    ax.tick_params(axis="x", labelsize=round(tick_fs * 0.8))
    for lbl in ax.get_yticklabels():
        lbl.set_fontweight(tick_fw)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    out_path = plots_dir / "gwas_phrog_functions.png"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  [phrog_bar] → {out_path.name}")
