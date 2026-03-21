# K-locus GRR — standalone runner

Computes all-vs-all weighted gene repertoire relatedness (wGRR) for 3,911 *Klebsiella pneumoniae* K-locus sequences. This is a long-running job (~hours); the computation is separated here so it can be run independently of the main chapter 2 analysis pipeline.

## Overview

1. Extracts CDS translations from each `.gb` file in `4_K_LOCI/` → per-isolate `.faa`
2. Merges all proteins into a single FASTA (proteins tagged `SAMPLE__PROTID`)
3. Builds one BLAST protein database
4. Runs a single all-vs-all `blastp` query
5. Identifies bidirectional best hits (BBH, e-value < 0.1) per pair
6. Aligns each BBH pair with `powerneedle` to get % identity
7. Computes wGRR = Σ(identity/100) / min(genes₁, genes₂)

## Input data

| Path | Description |
|------|-------------|
| `input_dir/gwas/4_K_LOCI/*.gb` | Per-isolate K-locus GenBank files (produced by `other/klocus_extraction`) |

## Output

| Path | Description |
|------|-------------|
| `output_dir/chapter2/grr_results.csv` | Pairwise wGRR table (columns: sample1, sample2, wgrr, grr_sum, genes1, genes2, bbh, min35, min80) |
| `output_dir/chapter2/grr_workdir/` | Intermediate files (faa/, BLAST DB, blast TSV); safe to delete once `grr_results.csv` exists |

## Estimated runtime (Apple Silicon, 8 threads)

| Step | Time |
|------|------|
| .faa extraction (3,911 isolates) | ~5 min |
| makeblastdb (~78 K proteins) | ~1 min |
| all-vs-all blastp | ~1–2 h |
| BBH parsing | ~5–10 min |
| powerneedle (pairs with ≥1 BBH) | ~2–6 h |
| **Total** | **~4–10 h** |

## Notes

- Requires `conda env jkoszucki` (BLAST+, Biopython) and `powerneedle` on PATH.
- Run via `conda run -n jkoszucki python klocus_grr.py`.
- Steps are skipped automatically if intermediate files already exist (resumable).
