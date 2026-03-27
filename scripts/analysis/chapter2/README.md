# Chapter 2 — GWAS analysis

Selects the best prophage-PC predictor per K-locus from pyseer GWAS results, exports
per-locus datasets (alignments, structures, GenBank files), and runs prophage BLAST.

---

## Pipeline steps (`main.py`)

| Step | Function | Output |
|------|----------|--------|
| 1 | `select_best_predictors` | `best_predictors.csv`, `best_predictors_functions.tsv` |
| 2 | `export_per_locus` | Alignment FASTAs + AF3 CIFs in `per-locus/{KL}/protein/` and `per-locus/experimental/{subfolder}/` |
| 3 | `blast_repr_vs_prophage` | Prophage BLAST results in `per-locus/{KL}/protein/` and `per-locus/experimental/{subfolder}/` |
| 4 | `export_per_protein` | Sequence FASTAs in `per-locus/{KL}/protein/` and `per-locus/experimental/{subfolder}/` |
| 5 | `export_klocus_genbank_dataset` | GenBank files sorted into `per-locus/{KL}/locus/` |
| 6 | `export_grr_vs_reference` | `per-locus/{KL}/locus/grr_vs_reference.tsv` (expensive — toggle with run flag) |

### Run flags (top of `main.py`)

```python
RUN_PROPHAGE_BLAST   = False   # True → runs blastp and overwrites raw TSV; False → skips blastp,
                                # still applies filter + MAFFT alignment on existing raw TSV
RUN_GRR_VS_REFERENCE = False   # True → runs BLAST + powerneedle per locus, overwrites grr_vs_reference.tsv
```

---

## Step 1 — Best predictor selection

### Filtering logic

1. From HHSEARCH annotation files, identify `(version, PC)` pairs with `tcov ≥ 0.10` and `reported_topology_PC == "SGNH hydrolase"`
2. Keep GWAS rows with `precision ≥ 0.6`
3. Inner-join with SGNH-qualified `(version, PC)` pairs
4. For each K-locus, retain the single `(PC, version)` with the highest `F1_score` across all clustering levels
5. Map representative sequence (first sequence of alignment FASTA)
6. Export all hhsearch rows for selected pairs (prob ≥ 0.70, no ALANDB)

### Input data

| File | Description |
|------|-------------|
| `input_dir/gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv` | All pyseer hits across clustering thresholds (~235k rows) |
| `input_dir/gwas/3_GWAS/1_INTERMEDIATE/3_FUNCTIONS/HHSEARCH/*/clusters_functions.tsv` | HHsearch annotations per clustering level (fold, tcov) |
| `input_dir/gwas/3_GWAS/1_INTERMEDIATE/2_MMSEQS/{version}/3_binary_matrix.tsv` | PC presence/absence matrix (rows=PCs, cols=genomeIDs) |

### Output

| File | Description |
|------|-------------|
| `output_dir/chapter2/best_predictors.csv` | One row per K-locus: locus, PC, version, mode, F1, representative_sequence, … |
| `output_dir/chapter2/best_predictors_functions.tsv` | HHsearch hits (prob ≥ 0.70) for selected predictors |

---

## Step 2 — Per-locus export (`export_per_locus`)

Copies alignment FASTAs and AF3 CIF structures into their final locations.

### Output structure

```
output_dir/chapter2/per-locus/
├── {KL}/
│   └── protein/
│       ├── {KL}_{PC}.fasta          — MMseqs2 alignment FASTA
│       └── {KL}_{PC}.cif            — AF3 model_0 structure
└── experimental/
    ├── 164_08_KL2_4NPA/
    │   └── 164_08_KL2.cif
    ├── 174_38_KL55_4NPA/
    │   └── 174_38_KL55.cif
    ├── QNO11465.1_K26/
    │   └── QNO11465.1_K26.cif
    └── YP_004782195.1_OAg_LPS/
        └── YP_004782195.1.cif
```

---

## Step 3 — Prophage BLAST (`blast_repr_vs_prophage`)

BLASTs each predictor's representative sequence and each experimental protein against all prophage proteins.

### Input

| Path | Description |
|------|-------------|
| `input_dir/gwas/2_PROPHAGES/4_FASTA_CDS_AA/*.faa` | Prophage protein FASTAs (~477k proteins) |
| `best_predictors.csv` | Representative sequences |
| `input_dir/gwas/sgnh-active.fasta` | Experimental protein sequences |

### Output (inside each protein's folder)

| File | Description |
|------|-------------|
| `{stem}_blast_raw.tsv` | All blastp hits (evalue ≤ 1e-5) |
| `{stem}_blast_filtered.tsv` | Hits with bitscore ≥ 500 |
| `{stem}_unaligned.fasta` | Query + filtered hit sequences |
| `{stem}_prophage.fasta` | MAFFT-aligned version (used by figures) |

Bitscore threshold: **500** (long aligned fragments only).

---

## Step 4 — Per-protein export (`export_per_protein`)

Writes sequence FASTAs for each predictor and experimental protein.

### Output

| File | Description |
|------|-------------|
| `per-locus/{KL}/protein/{KL}_{PC}_sequence.fasta` | Single-sequence FASTA for representative |
| `per-locus/experimental/{subfolder}/{stem}_sequence.fasta` | Single-sequence FASTA for experimental |

---

## Step 5 — K-locus GenBank dataset (`export_klocus_genbank_dataset`)

Classifies all Kpsc genomes per K-locus using two binary vectors:
- **locus** = K_locus typed as this type AND confidence ∈ {Perfect, Very high, High, Good}
- **PC** = best predictor PC present in genome (from binary matrix)

### Output structure

```
output_dir/chapter2/per-locus/{KL}/locus/
├── match/           ← locus==1 AND PC==1
├── locus_only/      ← locus==1 AND PC==0  (typed correctly, predictor absent)
├── pc_only/         ← locus==0 AND PC==1  (any Kpsc genome with PC present)
└── {KL}_reference.gb
```

All three categories span the **entire Kpsc dataset** (3900+ genomes). `pc_only` captures
genomes of any K-locus type in which the predictor fires but no good-confidence locus of
this type is detected.

### Input

| File | Description |
|------|-------------|
| `input_dir/gwas/1_BACTERIA/bacteria_metadata.tsv` | Kpsc genome metadata (K_locus, K_locus_confidence, species) |
| `input_dir/gwas/4_K_LOCI/{genomeID}.gb` | Per-genome K-locus GenBank files |
| `input_dir/cps/cps_loci_ref/K_locus_primary_referenece_MGG/{KL}.gb` | Reference GenBank per locus |
| `binary_matrix.tsv` (via mmseqs_dir) | PC presence/absence |

---

## Step 6 — wGRR vs reference (`export_grr_vs_reference`)

Computes wGRR similarity between the reference K-locus GenBank and every genome in
`match/`, `locus_only/`, and `pc_only/`. Uses GRRpair methodology: BBH via BLAST
(evalue < 0.1) + powerneedle + wGRR = Σ(identity/100) / min(genes_ref, genes_genome).

**Expensive** — runs BLAST + powerneedle for hundreds of genomes per locus. Toggle with
`RUN_GRR_VS_REFERENCE = True`.

### Output

| File | Description |
|------|-------------|
| `per-locus/{KL}/locus/grr_vs_reference.tsv` | genomeID, category, wgrr, grr_sum, genes_ref, genes_genome, n_bbh |

---

## Figures produced by `figures/chapter2/`

See `scripts/figures/chapter2/` for plot scripts. Key outputs:

| Plot | Script | Description |
|------|--------|-------------|
| `figure_sgnh_recall.png/.pdf` | `sgnh_recall_plot.py` | Recall per K-locus; marker shape = clustering level, size = PC abundance |
| `functions_heatmap.png` | `functions_heatmap.py` | PFAM domain presence/absence heatmap per best predictor (prob ≥ 0.90) |
| `gap_frequency_profile.png` | `alignment_profiles.py` | Gap-frequency profiles: MMseqs2 hits (top) vs prophage BLAST hits (bottom) |
| `per-locus/{KL}/{KL}.png` | `phandango_render.py` | Cladogram + K-locus confidence + PC presence heatmap (201-leaf subsampled tree) |
| `per-locus/{KL}/_variants.csv` | `phandango_export.py` | Phandango-ready CSV (filtered to tree leaves) |
| `per-locus/{KL}/{KL}_grr_vs_reference.png` | `grr_vs_reference_plot.py` | Box plot + jitter of wGRR to reference per category (match/locus_only/pc_only) |
| `per-locus/{KL}/{KL}_{PC}.png` | `structure_render.py` | AF3 CIF rendered to PNG via PyMOL (pymol conda env) |
| `per-locus/experimental/{subfolder}/{stem}.png` | `structure_render.py` | Experimental protein AF3 CIF rendered to PNG |
| `per-locus_legend.png` | `phandango_legend.py` | Colour legend for phandango PNGs |
| `sgnh_network_nodes.tsv` | `sgnh_network.py` | SGNH similarity network node table for Cytoscape |
| `sgnh_network_edges.tsv` | `sgnh_network.py` | SGNH similarity network edge table (all-vs-all BLAST) |
