"""
Gap-frequency profile plot for best-predictor protein alignments.

For each multi-FASTA (one per best predictor):
  1. Align with MAFFT (--auto, single-threaded, quiet).
  2. Compute per-column gap frequency.
  3. Bin into N_BINS equal-width bins along the normalised position (0–100%).
  4. Plot:
       - Individual labeled lines, one per alignment (K-locus + PC).
       - Median line and IQR ribbon across all alignments.

K-locus and PC are parsed from the filename convention  {KL-locus}_{PC}.fasta
(e.g. KL1_PC1817.fasta).

Run via:
    conda run -n alignment python alignment_profiles.py

Output:
    gap_frequency_profile.png
    gap_frequency_profile.svg
"""

from __future__ import annotations

import io
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Bio import SeqIO

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

_HELPERS_DIR = Path(__file__).resolve().parents[2] / "scripts" / "helpers"
sys.path.insert(0, str(_HELPERS_DIR))
from config import Config

cfg = Config()

FASTA_GLOB = cfg.output_dir / "chapter2" / "per-locus" / "*" / "*.fasta"
OUT_DIR    = Path(__file__).resolve().parent
N_BINS     = 30
GAP_CHARS  = set("-.")

# Palette: 18 distinct colours cycling through tab20
PALETTE = plt.cm.tab20.colors  # 20 colours


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_label(fasta_path: Path) -> str:
    """Return '{K-locus} {PC}' parsed from filename like KL1_PC1817.fasta."""
    stem = fasta_path.stem          # e.g. KL1_PC1817
    m = re.match(r"(KL\d+)_(PC\d+)", stem, re.IGNORECASE)
    if m:
        return f"{m.group(1)} {m.group(2)}"
    return stem


def _align(fasta_path: Path) -> list[str]:
    """Run MAFFT on a raw multi-FASTA; return list of aligned sequence strings."""
    result = subprocess.run(
        ["mafft", "--auto", "--quiet", str(fasta_path)],
        capture_output=True, text=True, check=True,
    )
    records = list(SeqIO.parse(io.StringIO(result.stdout), "fasta"))
    return [str(r.seq) for r in records]


def _gap_profile(seqs: list[str]) -> np.ndarray:
    """Return binned gap-frequency profile of shape (N_BINS,)."""
    aln_len = len(seqs[0])
    n_seqs  = len(seqs)
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
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    fasta_files = sorted(FASTA_GLOB.parent.parent.glob(
        str(FASTA_GLOB.relative_to(FASTA_GLOB.parent.parent))
    ))

    if not fasta_files:
        sys.exit(f"No FASTA files found matching {FASTA_GLOB}")

    profiles: list[np.ndarray] = []
    labels:   list[str]        = []

    for fasta in fasta_files:
        label = _parse_label(fasta)
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
            print(f"    [WARN] Unequal lengths in {fasta.name} after alignment, skipping.")
            continue
        profiles.append(_gap_profile(seqs))
        labels.append(label)
        print(f"    {len(seqs)} seqs, {lengths.pop()} cols → {label}")

    if not profiles:
        sys.exit("No profiles computed — nothing to plot.")

    mat        = np.array(profiles)           # (n, N_BINS)
    bin_edges  = np.linspace(0, 100, N_BINS + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    median     = np.nanmedian(mat, axis=0)
    q1         = np.nanpercentile(mat, 25, axis=0)
    q3         = np.nanpercentile(mat, 75, axis=0)
    n          = len(profiles)

    # --- Plot ---
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(10, 5))

    # Individual lines — pale gray, no legend entry
    for profile in profiles:
        ax.plot(bin_centers, profile, color="#cccccc", linewidth=0.8, alpha=0.7)

    # IQR ribbon — dark blue, clearly above gray
    ax.fill_between(bin_centers, q1, q3, color="#1f77b4", alpha=0.30,
                    label="IQR (25th–75th percentile)")

    # Median line — solid dark blue, thick
    ax.plot(bin_centers, median, color="#1f77b4", linewidth=2.5,
            linestyle="-", label="Median")

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Normalised alignment position (%)", fontsize=12)
    ax.set_ylabel("Gap frequency", fontsize=12)
    ax.set_title(f"Gap frequency profile across alignments (n={n})", fontsize=13)
    ax.tick_params(labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=10, frameon=True, loc="upper right")

    plt.tight_layout()
    for ext in ("png", "svg"):
        out = OUT_DIR / f"gap_frequency_profile.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else None, bbox_inches="tight")
        print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
