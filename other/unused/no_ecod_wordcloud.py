"""
Chapter 2, Figure 1, Panel B
─────────────────────────────
Word cloud of FoldSeek hit descriptions for no-ecod representative sequences.

Input:
    FoldSeek TSV results file (standard outfmt 0 / custom format with a
    'target_name' or 'description' column).  The file is placed manually by
    the user in {output_dir}/raw/ after running FoldSeek against the
    UniProt50/AlphaFold v6 database with taxonomic filter 'Bacteria'.

Expected columns (auto-detected):
    - 'target_name' or 'theader' — full description string for each hit
    If neither column is present, the second column is used as the description.

Frequency counting:
    Each unique description is split into words; common stopwords and short
    tokens are removed.  Word size in the cloud is proportional to frequency
    across all hits.

Outputs:
    plots_dir/figure1-panelB.png
    plots_dir/figure1-panelB.pdf
"""

from __future__ import annotations

import re
from pathlib import Path
from collections import Counter

import matplotlib.pyplot as plt
import pandas as pd


# Tokens to skip when counting words
_STOPWORDS = {
    "protein", "proteins", "domain", "domains", "family", "like",
    "homolog", "homologue", "related", "putative", "uncharacterized",
    "os", "ox", "gn", "pe", "sv",   # FASTA header tags
    "the", "a", "an", "of", "in", "and", "or", "for", "with",
    "to", "is", "are", "at", "from", "on", "by", "n",
}
_MIN_WORD_LEN = 3


def _parse_foldseek_tsv(foldseek_tsv: Path) -> pd.DataFrame:
    """
    Load a FoldSeek result TSV.  Handles both headered and headerless files.
    Returns a DataFrame with at least a 'description' column.
    """
    # Try reading with header first
    df = pd.read_csv(foldseek_tsv, sep="\t", low_memory=False)

    desc_col = None
    for candidate in ("target_name", "theader", "ttitle", "target"):
        if candidate in df.columns:
            desc_col = candidate
            break

    if desc_col is None:
        # Headerless file: assume standard FoldSeek outfmt 0
        # qseqid tseqid pident alnlen mismatch gapopen qstart qend tstart tend evalue bits
        df = pd.read_csv(foldseek_tsv, sep="\t", header=None)
        desc_col = 1  # target ID / description is typically the second column

    df = df.rename(columns={desc_col: "description"})
    return df[["description"]].dropna()


def _count_words(descriptions: pd.Series) -> Counter:
    """Tokenise descriptions and count word frequencies."""
    counts: Counter = Counter()
    for desc in descriptions:
        tokens = re.split(r"[\s\|_\-/;,]+", str(desc))
        for token in tokens:
            word = re.sub(r"[^a-zA-Z]", "", token).lower()
            if len(word) >= _MIN_WORD_LEN and word not in _STOPWORDS:
                counts[word] += 1
    return counts


def plot_no_ecod_wordcloud(
    foldseek_tsv: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Generate a word cloud from FoldSeek hits for no-ecod representative sequences.

    Args:
        foldseek_tsv:  path to FoldSeek TSV results file (placed by user in output_dir/raw/)
        plots_dir:     output directory for plots
        style:         cfg.style (optional)
    """
    try:
        from wordcloud import WordCloud
    except ImportError:
        print("  [figure1-panelB] wordcloud package not installed — skipping. "
              "Install with: pip install wordcloud")
        return

    if not foldseek_tsv.exists():
        print(f"  [figure1-panelB] FoldSeek results not found: {foldseek_tsv}\n"
              "  Place the file there manually and re-run with PLOT_FIGURE1_PANELB = True")
        return

    dpi = getattr(style, "dpi", 300)

    print(f"  [figure1-panelB] Loading FoldSeek results from {foldseek_tsv.name} …")
    df      = _parse_foldseek_tsv(foldseek_tsv)
    counts  = _count_words(df["description"])
    print(f"  [figure1-panelB] {len(df)} hits, {len(counts)} unique words after filtering")

    wc = WordCloud(
        width=1200,
        height=600,
        background_color="white",
        colormap="tab10",
        prefer_horizontal=0.9,
        min_font_size=8,
        max_words=150,
        collocations=False,
    ).generate_from_frequencies(counts)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.imshow(wc, interpolation="bilinear")
    ax.axis("off")
    plt.tight_layout(pad=0)

    for ext in ("png", "pdf"):
        out = plots_dir / f"figure1-panelB.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"  [figure1-panelB] → {out.name}")
    plt.close(fig)
