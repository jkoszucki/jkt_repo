"""
Gap-frequency profile plot for best-predictor protein alignments.

Two subplots sharing the X-axis are produced when `prophage_blast_dir` is
supplied; otherwise a single subplot is produced (backward-compatible).

Top subplot (always):
  For each best predictor (locus, PC):
    1. Locate the multi-FASTA: per_locus_dir / {locus} / protein / {locus}_{PC}.fasta
    2. Align with MAFFT (--auto, quiet, in-memory).
    3. Compute per-column gap frequency; bin into N_BINS equal-width bins
       over normalised position 0–100%.
  Plot individual pale-gray lines + median + IQR ribbon.

Bottom subplot (when prophage_blast_dir is given):
  Reads pre-computed aligned FASTAs from
      prophage_blast_dir / {KL}_{PC}_prophage.fasta
  produced by analysis/chapter2 blast_repr_vs_prophage().
  Computes gap profiles and shows in the same way.

Font sizes and weights follow cfg.style (same as functions_heatmap).
"""

from __future__ import annotations

import io
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO

N_BINS    = 30
GAP_CHARS = set("-.")


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def _align(fasta_path: Path) -> list[str]:
    result = subprocess.run(
        ["mafft", "--auto", "--quiet", str(fasta_path)],
        capture_output=True, text=True, check=True,
    )
    return [str(r.seq) for r in SeqIO.parse(io.StringIO(result.stdout), "fasta")]


def _gap_profile(seqs: list[str]) -> np.ndarray:
    aln_len  = len(seqs[0])
    n_seqs   = len(seqs)
    gap_freq = np.array([
        sum(1 for s in seqs if s[i] in GAP_CHARS) / n_seqs
        for i in range(aln_len)
    ])
    positions  = np.linspace(0, 100, aln_len)
    bin_edges  = np.linspace(0, 100, N_BINS + 1)
    bin_means  = []
    for b in range(N_BINS):
        lo, hi = bin_edges[b], bin_edges[b + 1]
        mask = (positions >= lo) & (positions <= hi if b == N_BINS - 1 else positions < hi)
        bin_means.append(gap_freq[mask].mean() if mask.sum() > 0 else np.nan)
    return np.array(bin_means)


# ---------------------------------------------------------------------------
# Profile computation
# ---------------------------------------------------------------------------

def _compute_profiles_from_fastas(
    best: pd.DataFrame, per_locus_dir: Path
) -> list[np.ndarray]:
    """Align per-locus MMseqs2 FASTAs and compute gap profiles."""
    profiles: list[np.ndarray] = []
    for _, row in best.iterrows():
        locus, pc = row["locus"], row["PC"]
        fasta = per_locus_dir / locus / "protein" / f"{locus}_{pc}.fasta"
        if not fasta.exists():
            print(f"  [alignment] missing {fasta.name}, skipping.")
            continue
        print(f"  Aligning {fasta.name} …", flush=True)
        try:
            seqs = _align(fasta)
        except subprocess.CalledProcessError as e:
            print(f"    [WARN] MAFFT failed for {fasta.name}: {e}")
            continue
        if len(seqs) < 2:
            print(f"    [WARN] < 2 sequences in {fasta.name}, skipping.")
            continue
        lengths = {len(s) for s in seqs}
        if len(lengths) > 1:
            print(f"    [WARN] unequal lengths after alignment in {fasta.name}, skipping.")
            continue
        profiles.append(_gap_profile(seqs))
        print(f"    {len(seqs)} seqs, {lengths.pop()} cols → {locus} {pc}")
    return profiles


def _compute_profiles_from_prophage_blast(
    best: pd.DataFrame, per_locus_dir: Path
) -> list[np.ndarray]:
    """Read pre-computed prophage-aligned FASTAs and compute gap profiles."""
    profiles: list[np.ndarray] = []
    for _, row in best.iterrows():
        locus, pc = row["locus"], row["PC"]
        fasta = per_locus_dir / locus / "protein" / f"{locus}_{pc}_prophage.fasta"
        if not fasta.exists():
            print(f"  [prophage] missing {fasta.name}, skipping.")
            continue
        seqs = [str(r.seq) for r in SeqIO.parse(fasta, "fasta")]
        if len(seqs) < 2:
            print(f"  [prophage] < 2 sequences in {fasta.name}, skipping.")
            continue
        lengths = {len(s) for s in seqs}
        if len(lengths) > 1:
            print(f"  [prophage] unequal lengths in {fasta.name}, skipping.")
            continue
        profiles.append(_gap_profile(seqs))
        print(f"  {fasta.name}: {len(seqs)} seqs, {lengths.pop()} cols")
    return profiles


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _draw_profile_ax(
    ax,
    profiles: list[np.ndarray],
    title: str,
    label_fs: int,
    label_fw: str,
    tick_fs: int,
    tick_fw: str,
    color: str = "#1f77b4",
) -> None:
    bin_edges   = np.linspace(0, 100, N_BINS + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    mat         = np.array(profiles)
    median      = np.nanmedian(mat, axis=0)
    q1          = np.nanpercentile(mat, 25, axis=0)
    q3          = np.nanpercentile(mat, 75, axis=0)
    n           = len(profiles)

    for profile in profiles:
        ax.plot(bin_centers, profile, color="#cccccc", linewidth=0.8, alpha=0.30)

    ax.fill_between(bin_centers, q1, q3, color=color, alpha=0.30,
                    label="IQR (25th–75th percentile)")
    ax.plot(bin_centers, median, color=color, linewidth=2.5,
            linestyle="-", label="Median")

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("Gap frequency", fontsize=label_fs, fontweight=label_fw)
    ax.set_title(f"{title} (n={n}, bins={N_BINS})", fontsize=label_fs + 1, fontweight=label_fw)
    ax.tick_params(labelsize=tick_fs)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight(tick_fw)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=tick_fs, frameon=True, loc="upper right")


def plot_alignment_profiles(
    per_locus_dir: Path,
    best_predictors_csv: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Align best-predictor FASTAs and plot gap-frequency profiles.

    Args:
        per_locus_dir:       cfg.output_dir / "chapter2" / "per-locus"
        best_predictors_csv: cfg.output_dir / "chapter2" / "best_predictors.csv"
        plots_dir:           scripts/figures/chapter2/plots
        style:               cfg.style
    """
    tick_fs  = getattr(style, "tick_fontsize",        8)
    tick_fw  = getattr(style, "tick_fontweight",      "bold")
    label_fs = getattr(style, "axis_label_fontsize",  12)
    label_fw = getattr(style, "axis_label_fontweight","bold")
    dpi      = getattr(style, "dpi",                  300)

    best = pd.read_csv(best_predictors_csv)

    top_profiles = _compute_profiles_from_fastas(best, per_locus_dir)
    if not top_profiles:
        print("  [alignment] no profiles computed — skipping plot.")
        return

    bot_profiles = _compute_profiles_from_prophage_blast(best, per_locus_dir)

    sns.set_style("whitegrid")
    two_panels = bool(bot_profiles)

    if two_panels:
        fig, (ax_top, ax_bot) = plt.subplots(
            2, 1, figsize=(10, 9), sharex=True,
            gridspec_kw={"hspace": 0.35},
        )
    else:
        fig, ax_top = plt.subplots(figsize=(10, 5))
        ax_bot = None

    _draw_profile_ax(
        ax_top, top_profiles,
        title="MMseqs2 hit alignments",
        label_fs=label_fs, label_fw=label_fw,
        tick_fs=tick_fs, tick_fw=tick_fw,
        color="#1f77b4",
    )

    if two_panels:
        _draw_profile_ax(
            ax_bot, bot_profiles,
            title="Prophage BLAST hit alignments",
            label_fs=label_fs, label_fw=label_fw,
            tick_fs=tick_fs, tick_fw=tick_fw,
            color="#d62728",
        )
        ax_bot.set_xlabel(
            "Normalised alignment position (%)",
            fontsize=label_fs, fontweight=label_fw,
        )
    else:
        ax_top.set_xlabel(
            "Normalised alignment position (%)",
            fontsize=label_fs, fontweight=label_fw,
        )

    plt.tight_layout()
    out_path = plots_dir / "gap_frequency_profile.png"
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out_path}")
