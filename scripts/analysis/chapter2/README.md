# Chapter 2 — GWAS analysis

Selects the best prophage-PC predictor per K-locus from pyseer GWAS results.

## Filtering logic

1. Keep rows with `precision ≥ 0.8`
2. For each K-locus, retain the single `(PC, version)` combination with the highest `F1 × MCC` product across all clustering thresholds (`version`)
3. Keep only loci where `F1 ≥ 0.5` and `MCC ≥ 0.5`

## Input data

| File | Description |
|------|-------------|
| `input_dir/gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv` | All pyseer hits across clustering thresholds (~235k rows) |

## Output

| File | Description |
|------|-------------|
| `output_dir/chapter2/best_predictors.csv` | One row per K-locus passing all filters |
| `output_dir/chapter2/grr_results.csv` | All-vs-all pairwise wGRR for 3,911 K-loci (columns: sample1, sample2, wgrr, grr_sum, genes1, genes2, bbh, min35, min80) |
| `output_dir/chapter2/grr_workdir/` | Intermediate files (per-isolate .faa, merged FASTA, BLAST DB, BLAST TSV); safe to delete after `grr_results.csv` is written |
