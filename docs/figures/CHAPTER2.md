# Chapter 2 — GWAS: proteins without annotated folds are atypical tail fiber proteins and acetyltransferases

Chapter 2 focuses on GWAS-predicted prophage proteins that **lack an assigned ECOD fold** (no-ecod proteins), showing they are atypical tail fiber proteins or acetyltransferases. It also covers prophage-encoded acetyltransferases detected by HHsearch and their co-occurrence with K-loci.

> **Data preparation** is split across two processing modules — see `docs/processing/GWAS-PROC.md` and `docs/processing/ACETYL-PROC.md` for pipeline detail, inputs, outputs, and library modules. `figures/chapter2/` is plotting only: it reads from `output_dir/gwas-data/`, `output_dir/acetyltransferase/`, and `output_dir/manual-outputs/` and writes plots to `scripts/figures/chapter2/plots/`.

---

## Chapter 2 specific scope

| Analysis | Handled by | Description |
|----------|-----------|-------------|
| No-ecod representative clustering | `figures/chapter2/` | Cluster no-ecod representative sequences with MMseqs2 (e-value ≤ 10⁻³, id ≥ 80%, cov ≥ 80%); select one representative per cluster |
| Foldseek search (no-ecod) | manual → `output_dir/manual-outputs/foldseek/` | Predict AF3 monomers for no-ecod representatives; search vs UniProt50/AlphaFold v6 (Bacteria filter) |
| No-ecod similarity network | `figures/chapter2/` | All-vs-all BLASTP (e-value ≤ 10⁻³, id ≥ 50%, cov ≥ 80%); edge.tsv + node.tsv for Cytoscape |
| SSLBH GWAS predictor detection | `processing/acetyl-proc/` | From `acetyl-gwas/pyseer_hits_sslbh_pvalcor005.tsv`: filter `sslbh_detected == True`; select best predictor per K-locus (highest F1, precision ≥ 0.60) |
| Prophage acetyltransferase detection | `processing/acetyl-proc/` | HHsearch at PC80 level (ECOD SSLBH prob ≥ 0.7; PHROGs acetyltransferase keywords prob ≥ 0.9; PFAM PF01757 prob ≥ 0.9); all detected PC80s in `acetyl-annot-pc80/` |
| Experimental acetyltransferases | manual → `output_dir/manual-outputs/` | Foldseek vs UniProt50/AF v6; HHphred vs PFAM and ECOD for 3 experimental proteins |
| TM-align screen (no-ecod GWAS) | `processing/acetyl-proc/` | TM-align AF3 monomers of no-ecod GWAS predictors vs reference acetyltransferases; results in `acetyltransferase/tmalign_results.tsv` |

---

## Input data

Data prepared by `processing/gwas-proc/`:

| File | Description |
|------|-------------|
| `output_dir/gwas-data/gwas_hits.tsv` | All GWAS hits passing quality filters — filter for `ecod_type == "no-ecod"` for tail fiber / TM-align analysis |
| `output_dir/gwas-data/{ecod-type}/{k-locus}/{clustering-level}/{PC}/protein/` | Protein data per predictor (sequence.fasta, structure.cif, raw_blast.tsv) |

Data prepared by `processing/acetyl-proc/`:

| File | Description |
|------|-------------|
| `output_dir/acetyltransferase/acetyl-gwas/pyseer_hits_sslbh_pvalcor005.tsv` | All pyseer hits (pvalue_corr ≤ 0.05) with `sslbh_detected`, `sslbh_pc80`, `sslbh_detected_by` columns — source for Figure 2B SSLBH predictor selection |
| `output_dir/acetyltransferase/acetyl-annot-pc80/acetyl-pc80.tsv` | All prophage-encoded acetyltransferase PC80s with protein → prophage → genome mapping — source for Figure 3A/3B |
| `output_dir/acetyltransferase/acetyl-annot-pc80/acetyl-hhsearch-pc80.tsv` | HHsearch hits that triggered PC80 detection — used for evidence annotation |
| `output_dir/acetyltransferase/acetyl-annot-pc80/{pc80}/` | Per-PC80: sequence.fasta, sequence.cif (if available), against-prophages/raw_blast.tsv |
| `output_dir/acetyltransferase/tmalign_results.tsv` | TM-align scores for no-ecod GWAS predictors vs reference acetyltransferases |

Chapter 2 specific inputs:

| File | Description |
|------|-------------|
| `input_dir/cps/cps.xlsx` (sheet: `modifying_enzymes`) | Experimental acetyltransferase metadata |
| `output_dir/manual-outputs/alphafold3/` | AF3 predictions (placed manually) |
| `output_dir/manual-outputs/foldseek/` | Foldseek results (placed manually) |

Experimental acetyltransferases: `PROTEIN01_MOD_AC_K1`, `PROTEIN02_MOD_AC_K2`, `PROTEIN03_MOD_AC_K57`

---

## Figures produced by `figures/chapter2/`

| Figure | Panel | Plot file | Description |
|--------|-------|-----------|-------------|
| Figure 1 | A | `figure1-panelA.png` | 6-panel scatter (F1 vs MCC) per clustering level, coloured by ECOD topology |
| Figure 1 | B | `figure1-panelB.png` | Word cloud from FoldSeek hits (no-ecod representatives vs UniProt50/AlphaFold v6) |
| Figure 1 | C | `figure1-panelC.png` | Sequence similarity network of no-ecod proteins (Cytoscape edge.tsv + node.tsv) |
| Figure 2 | A | `figure2-panelA.png` | Experimental acetyltransferases: FoldSeek word cloud + ECOD/PFAM annotation |
| Figure 2 | B | `figure2-panelB.png` | GWAS SSLBH predictors per K-locus (dot/bar plot) |
| Figure 3 | A | `figure3-panelA.png` | Stacked bar: acetyltransferase distribution across prophages per K-locus |
| Figure 3 | B | `figure3-panelB.png` | Diversity network: prophage-encoded + GWAS + experimental acetyltransferases (Cytoscape) |
| Figure 4 | A | `figure4-panelA.png` | Clinker-style gene maps of prophage acetyltransferase genomic context |

---

## Methodology notes

**Figure 1A — ECOD topology colours:**
Each PC coloured by `ecod_type` and `reported_topology_PC` from `gwas_hits.tsv` (quality-filtered hits from `processing/gwas-proc`). Colour mapping:
- SGNH hydrolase → `sgnh_domain_color`
- SSRBH / Pectin-lyase like → `ssrbh_color`
- SSLBH → `sslbh_color`
- No ECOD hit → `no_ecod_color` (white fill, gray border)
- Other ECOD → `gray_color`

**Figure 2B — SSLBH GWAS predictor selection:**
Source: `acetyltransferase/acetyl-gwas/pyseer_hits_sslbh_pvalcor005.tsv`, filter `sslbh_detected == True`.
Best predictor per K-locus = highest F1 among rows with `sslbh_detected == True` and `precision ≥ 0.60`, across all clustering levels.

**Figure 3A — Prophage acetyltransferase detection:**
Source: `acetyltransferase/acetyl-annot-pc80/acetyl-pc80.tsv`.
Detection method: HHsearch at PC80 level with three evidence sources (ECOD SSLBH, PHROGs keywords, PFAM PF01757) — see `docs/processing/ACETYL-PROC.md` for thresholds. All 133 detected PC80s are included. K-locus assignment per PC80: if ≥ 50% of member sequences come from genomes with the same K-locus → assign that locus; otherwise `MANY-LOCI(K2,K13)`.

**Figure 1B/C — No-ecod clustering:**
Cluster representative sequences from no-ecod PCs using MMseqs2 easy-cluster (e-value ≤ 10⁻³, id ≥ 80%, cov ≥ 80%). Select one representative per cluster for AF3 prediction and FoldSeek search.
