# Alignment gap-frequency profiles

Standalone script that aligns multi-FASTA files with MAFFT and visualises per-alignment gap-frequency profiles for the 18 best GWAS predictors.

## What it does

1. Discovers all `*.fasta` files under `cfg.output_dir/chapter2/per-locus/*/`.
2. Aligns each multi-FASTA with `mafft --auto --quiet` (in-memory; no files written).
3. Computes per-column gap frequency for each alignment.
4. Bins each profile into 50 equal-width bins over normalised position 0–100%.
5. Plots:
   - One labelled line per alignment (`{K-locus} {PC}`, e.g. `KL1 PC1817`).
   - Dashed black median line across all alignments.
   - Gray IQR (25th–75th percentile) ribbon.

## Input

| Path | Description |
|------|-------------|
| `output_dir/chapter2/per-locus/{KL}/*.fasta` | Raw multi-FASTA files (one per best predictor, produced by `analysis/chapter2`) |

## Output

| File | Description |
|------|-------------|
| `gap_frequency_profile.png` | 300 dpi raster figure |
| `gap_frequency_profile.svg` | Vector figure |

Both files are written next to this script.

## Environment

Requires a dedicated `alignment` conda environment (MAFFT is not in `jkoszucki`):

```bash
conda create -n alignment -c conda-forge -c bioconda \
    mafft biopython matplotlib numpy seaborn pyyaml -y
```

Run:
```bash
cd other/alignment_profiles
conda run -n alignment python alignment_profiles.py
```
