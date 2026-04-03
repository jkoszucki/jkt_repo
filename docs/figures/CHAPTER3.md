# Chapter 3 — GWAS: SGNH hydrolase proteins are putative capsular polysaccharide deacetylases

Chapter 3 focuses on GWAS-predicted prophage proteins carrying the **SGNH hydrolase** ECOD fold, showing they are putative capsular polysaccharide deacetylases that co-occur with K-loci carrying documented O-acetylation.

> **Data preparation** (filtering GWAS hits, ECOD classification, per-locus protein/locus/prophage export) is done by `processing/gwas-proc/`. Chapter 3 reads from `output_dir/gwas-data/` filtered for `ecod_type == "sgnh-ecod"`.

---

## Chapter 3 specific scope

| Analysis | Description |
|----------|-------------|
| SGNH similarity network | All-vs-all BLASTP of 18 SGNH predictors + 4 experimental SGNH proteins (cov ≥ 30%, id ≥ 50%, e-value ≤ 10⁻³); edge.tsv + node.tsv for Cytoscape |
| Phandango export | Phylogenetic subtree + K-locus confidence + PC presence/absence per K-locus; exported by `processing/gwas-proc/` to `gwas-data/sgnh-ecod-reported-topology/{k-locus}/{clustering-level}/{PC}/phandango/`; `figures/chapter3/` reads these files and renders PNGs |
| Alignment profiles | Gap-frequency profiles from MMseqs2 FASTA and BLASTP-vs-prophages FASTA per predictor |

---

## Input data

Data prepared by `processing/gwas-proc/`:

| File | Description |
|------|-------------|
| `output_dir/gwas-data/gwas_sgnh_hits.tsv` | Best SGNH predictor per K-locus |
| `output_dir/gwas-data/sgnh-ecod-reported-topology/{k-locus}/{clustering-level}/{PC}/protein/` | Protein data per predictor |
| `output_dir/gwas-data/sgnh-ecod-reported-topology/{k-locus}/{clustering-level}/{PC}/locus/` | K-locus GenBanks and grr.tsv per predictor |

Chapter 3 specific inputs:

| File | Description |
|------|-------------|
| `input_dir/gwas/sgnh-active.fasta` | 4 experimentally characterised SGNH hydrolase sequences |
| `input_dir/gwas/1_BACTERIA/bacteria_metadata.tsv` | Kpsc genome metadata (for phandango phylogenetic subtree) |
| `output_dir/manual-outputs/alphafold3/` | AF3 trimer predictions (placed manually) |
| `output_dir/cps-structures/ktypes_modifications.csv` | Modification data for Figure 4A (sourced from chapter 4 processing) |

---

## Outputs

Written to `output_dir/gwas-data/sgnh-ecod-reported-topology/{k-locus}/{clustering-level}/{PC}/` by `processing/gwas-proc/`:

| Subpath | Description |
|---------|-------------|
| `phandango/subtree.nwk` | Phylogenetic subtree for this K-locus |
| `phandango/variants.csv` | PC presence/absence per genome for phandango rendering |

Experimental SGNH proteins written to `output_dir/gwas-data/experimental/{proteinID}/` by `processing/gwas-proc/`:

| Subpath | Description |
|---------|-------------|
| `structure.cif` | AF3 trimer (model_0) |
| `sequence.fasta` | Protein sequence |
| `against-prophages/raw_blast.tsv` | Unfiltered BLASTP hits vs prophage proteins |

---

## Figures produced by `figures/chapter3/`

| Figure | Panel | Plot file | Description |
|--------|-------|-----------|-------------|
| Figure 1 | A | `figure1-panelA.png` | Horizontal dot/bar plot: F1 score per K-locus for SGNH predictors |
| Figure 1 | B | `figure1-panelB.png` | Sequence similarity network of SGNH predictors + experimental proteins (Cytoscape edge.tsv + node.tsv) |
| Figure 2 | A | `figure2-panelA/` | Phandango PNGs per K-locus (phylogenetic subtree + K-locus + PC presence) |
| Figure 2 | B | `figure2-panelB.png` | GRR distribution histograms per category (match / locus_only / pc_only) across all K-loci |
| Figure 3 | A | `figure3-panelA.png` | PFAM domain profile rectangles per predictor sequence |
| Figure 3 | B | `figure3-panelB.png` | Gap-frequency profiles (median ± IQR, 30 bins) + conservation score column |
| Figure 4 | A | `figure4-panelA.png` | CPS repeating unit structures for K2, K8, K11, K24, K55 with O-acetylation annotations |
| Figure 4 | B | `figure4-panelB.png` | Clinker-style gene map of K8 prophage encoding the SGNH hydrolase predictor |

---

## Methodology notes

**Figure 1A — SGNH predictor selection:**
Best predictor per K-locus = PC with precision ≥ 60%, SGNH hydrolase ECOD topology (`reported_topology_PC == "SGNH hydrolase"`), tcov ≥ 10%, prob ≥ 70%, highest F1 across all six clustering levels. Representative protein sequence = first sequence from the MMseqs2 alignment FASTA. Reuse and update code from `figure_sgnh_recall.png`.

**Figure 1A — Visualisation:**
Horizontal dot/bar plot; x-axis = F1 score; rows = K-locus labels; colour = `sgnh_domain_color`; dot outline thickness distinguishes predictor quality: good (recall ≥ 0.6) → thick; likely (recall < 0.5) → thin. Sorted by K-locus number.

**Figure 1B — Sequence similarity network:**
All-vs-all BLASTP of 18 SGNH predictor representative sequences + 4 experimental SGNH hydrolases. Retain hits: cov ≥ 30%, id ≥ 50%, e-value ≤ 10⁻³. Coverage interval label per edge: low (30–50%), medium (50–80%), high (>80%). AF3 trimer structures placed into network manually. Reuse code from `sgnh_network_edges.tsv` / `sgnh_network_nodes.tsv` generation.

**Figure 2A — Phandango:**
Retain all K-locus confidence levels (not just "Good") for the phylogenetic subtree. One PNG per K-locus. Reuse and update code from phandango plots in `per-locus/`.

**Figure 2B — GRR distributions:**
Three-subplot column histogram (match / locus_only / pc_only), wGRR values from `locus/grr.tsv`, pooled across all 18 K-loci. Reuse and update code from `grr_vs_reference_global.png`.

**Figure 3A — PFAM domain profiles:**
Reuse snippet from `other/arch/draw_alignment.py` (originally used for a different purpose).

**Figure 3B — Gap-frequency and conservation:**
Two alignment sources per predictor: (1) MMseqs2 PC FASTA (`protein/pc.fasta`), (2) BLASTP hits vs prophages (`protein/against-prophages/`). Bins = 30; overlay = median ± IQR across all 18 predictors (no individual lines). Conservation score as subplot in second column. Reuse and update code from `gap_frequency_profile.png`.

**Figure 4A — CPS structures:**
Structures for K2, K8, K11, K24, K55 reused from `output_dir/cps-structures/` (generated by `processing/cps-proc/`). Annotate O-acetylation groups and occupancy using `ktypes_modifications.csv`.
