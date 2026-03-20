# Chapter 2 — K-locus extraction

Extracts the K-locus genomic region from annotated GenBank assemblies for 3,911 *Klebsiella pneumoniae* genomes.

## Overview

**Phase 1** — Kaptive v3 (`kaptive assembly kpsc_k`) is run on all FASTA assemblies to identify K-locus coordinates (contig, start, end).

**Phase 2** — For each sample, the detected region is sliced from the corresponding GenBank file using Biopython, preserving all CDS annotations within the locus. Multi-contig loci are concatenated. Samples with no detected locus produce a 0-byte placeholder file.

## Input data

| Path | Description |
|------|-------------|
| `input_dir/gwas/1_BACTERIA/1_GENBANK_GENOMES/*.gb` | Prokka-annotated GenBank files (one per assembly) |
| `input_dir/gwas/1_BACTERIA/2_FASTA_GENOMES_NT/*.fasta` | FASTA assemblies (one per assembly) |

## Output

| Path | Description |
|------|-------------|
| `output_dir/chapter2/kaptive_results.jsonl` | Kaptive JSON Lines output (one record per sample) |
| `input_dir/gwas/4_K_LOCI/*.gb` | Extracted K-locus GenBank files; 0-byte if no locus detected |

The LOCUS name in each output `.gb` file encodes the sample and K-locus type (e.g., `1001KBV_KL22`).

## Notes

- Requires `conda env kaptive` (Kaptive v3). Run via `conda run -n kaptive python main.py`.
- K-loci output is written to `input_dir/gwas/4_K_LOCI/` rather than `output_dir` because it serves as input for downstream GWAS analysis.
