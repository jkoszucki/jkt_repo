# processing/gwas-proc — GWAS data preparation

Shared data preparation module feeding **chapters 2 and 3**. No chapter-specific logic lives here — all outputs are consumed by `figures/chapter2/` and `figures/chapter3/` (plotting only).

## Pipeline

```
pyseer_hits_all.tsv  +  clusters_functions_best_all.tsv (provides reported_topology_PC per clustering_level × PC)
        │
        ▼ Step 1: filter by quality thresholds; join and classify by reported_topology_PC
gwas_hits.tsv   (all hits passing thresholds, one row per locus × ecod-type × clustering-level × PC)
        │
        ├─▶ Step 2: per-PC file export
        │   2a/2b → protein/pc.fasta, sequence.fasta
        │   2c    → output_dir/manual-upload/gwas_batch_NNN.json  (AF3 upload)
        │   2d/2e → prophage DB + against-prophages/raw_blast.tsv
        │   2f/2g → locus/reference.gb, grr.tsv, match/, locus_only/, pc_only/
        │
        └─▶ Step 3: select best SGNH predictors per K-locus
            → gwas_sgnh_hits.tsv (highest F1 from sgnh-ecod-reported-topology/, precision ≥ 0.60)
```

## Run flags

Defined at the top of `main.py`:

| Flag | Default | Controls |
|------|---------|----------|
| `RUN_AF3_JSON_PREP` | `True` | Write `gwas_batch_NNN.json` upload files; skips PCs where `structure.cif` already exists |
| `RUN_PROPHAGE_BLAST` | `True` | BLASTP of representative sequences vs prophage proteins; skips per-PC if `raw_blast.tsv` exists |
| `RUN_GRR_VS_REFERENCE` | `True` | wGRR computation (BLAST + powerneedle, expensive); skips per-PC if `grr.tsv` exists |

## Input data

| Path | Description |
|------|-------------|
| `input_dir/gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv` | All pyseer GWAS hits |
| `input_dir/gwas/3_GWAS/1_INTERMEDIATE/3_FUNCTIONS/clusters_functions_best_all.tsv` | Best HHsearch annotation per (PC, clustering_level); provides `reported_topology_PC` |
| `input_dir/gwas/3_GWAS/1_INTERMEDIATE/2_MMSEQS/{clustering_level}/` | MMseqs2 clustering outputs per level |
| `input_dir/gwas/2_PROPHAGES/4_FASTA_CDS_AA/*.faa` | Prophage protein FASTAs (for BLAST DB) |
| `input_dir/gwas/1_BACTERIA/bacteria_metadata.tsv` | Kpsc genome metadata |
| `input_dir/cps/cps_loci_ref/K_locus_primary_referenece_MGG/` | Reference K-locus GenBanks |
| `input_dir/gwas/4_K_LOCI/` | Per-genome K-locus GenBanks |
| `output_dir/manual-outputs/alphafold3/` | AF3 predictions (placed manually per batch) |

## Output data structure

Written to `output_dir/gwas-data/` and `output_dir/manual-upload/`. The prophage BLAST DB is shared and built separately by `scripts/helpers/prophage_db.py` (see ENZYMES-PROC.md).

```
manual-upload/
    gwas_batch_001.json              — AF3 upload batch (≤30 predictions per file)
    gwas_batch_002.json
    ...
    af3_manifest.tsv                 — af3_id | run_as | batch | sequence
    af3_placement.tsv                — af3_id | source_path | locus | ecod_type | clustering_level | PC | protein_id

gwas-data/
    gwas_hits.tsv                        — all hits passing thresholds; one row per (locus, ecod-type, clustering-level, PC)
    gwas_sgnh_hits.tsv                   — best SGNH predictor per K-locus (highest F1, precision ≥ 0.60)

    sgnh-ecod-reported-topology/{k-locus}/{clustering-level}/{PC}/
    ssrbh-ecod-reported-topology/{k-locus}/{clustering-level}/{PC}/
    other-ecod-reported-topology/{k-locus}/{clustering-level}/{PC}/
    no-ecod-reported-topology/{k-locus}/{clustering-level}/{PC}/
        protein/
            pc.fasta                     — prophage protein cluster as multifasta (NOT alignment)
            sequence.fasta               — representative sequence (first from alignment FASTA)
            structure.cif                — AF3 trimer (sgnh/ssrbh) or monomer (other/no-ecod); model_0
            against-prophages/
                raw_blast.tsv            — BLASTP hits vs prophage DB; sseqid = protein_name extracted
                                           as parts[1] from '|||'-delimited header ('>CDS_1 ||| protein_name ||| ...')
                                           columns: qseqid, sseqid, pident, length, mismatch, gapopen,
                                                    qstart, qend, sstart, send, evalue, bitscore,
                                                    qcovhsp, qlen, slen, sseq
        locus/
            reference.gb                 — reference K-locus GenBank from Kaptive database
            grr.tsv                      — wGRR vs reference; columns: genomeID, category, wgrr, ...
            match/                       — K-locus GenBanks: K-locus and PC co-occur
            locus_only/                  — K-locus GenBanks: K-locus present, PC absent
            pc_only/                     — K-locus GenBanks: PC present, K-locus not typed

```

## ECOD types

Classification uses `reported_topology_PC` joined from `clusters_functions_best_all.tsv` by (clustering_level, PC). Multiple rows per pair always carry the same topology value; one representative row is taken per pair.

| Folder | Topology | Classification rule |
|--------|----------|---------------------|
| `sgnh-ecod-reported-topology` | SGNH hydrolase | "sgnh" in `reported_topology_PC` |
| `ssrbh-ecod-reported-topology` | SSRBH — depolymerases | "pectin" in `reported_topology_PC` ("Pectin-lyase like" is the ECOD name for the SSRBH fold) |
| `other-ecod-reported-topology` | Other ECOD topology | Any annotated topology not matched above |
| `no-ecod-reported-topology` | No ECOD hit | `reported_topology_PC` absent, empty, or `"no hit"` |

## Quality filters (Step 1)

| Threshold | Value |
|-----------|-------|
| mode | `lasso` |
| pvalue_corr | ≤ 0.05 |
| Precision | ≥ 0.60 |
| F1 score | ≥ 0.50 |
| MCC | ≥ 0.50 |

All clustering levels are searched; every (PC, clustering_level) combination passing all thresholds is retained in `gwas_hits.tsv`.

Best SGNH predictor per K-locus = highest F1 from `sgnh-ecod-reported-topology/`, precision ≥ 0.60, across all clustering levels.

## Clustering levels

| Level | Identity | Coverage |
|-------|----------|----------|
| `PCI00C50` | 0% | 50% |
| `PCI50C50` | 50% | 50% |
| `PCI80C50` | 80% | 50% |
| `PCI00C80` | 0% | 80% |
| `PCI50C80` | 50% | 80% |
| `PCI80C80` | 80% | 80% |

## Library modules (`lib/`)

| Module | Function | Role |
|--------|----------|------|
| `gwas_filter_classify.py` | `filter_and_classify()` | Step 1: filter by thresholds, classify by `reported_topology_PC`; writes `gwas_hits.tsv` |
| `per_locus_export.py` | `export_per_locus()` | Step 2a: copy alignment FASTAs (→ `pc.fasta`) and AF3 CIFs (→ `structure.cif`) |
| `per_locus_export.py` | `export_per_protein()` | Step 2b: write `sequence.fasta` per predictor |
| `helpers/af3_json_prep.py` | `prepare_af3_batches()` | Step 2c: write `gwas_batch_NNN.json` for pending sequences (checkpoint: skip if `structure.cif` exists) |
| `helpers/prophage_db.py` | `prepare_prophage_db()` | Step 2d: build shared BLAST DB from prophage FASTAs |
| `per_locus_export.py` | `blast_repr_vs_prophage()` | Step 2e: BLASTP vs prophage DB; writes `against-prophages/raw_blast.tsv` |
| `per_locus_export.py` | `export_klocus_genbank_dataset()` | Step 2f: copy reference GenBank + classify genomes into match/locus_only/pc_only |
| `per_locus_grr.py` | `export_grr_vs_reference()` | Step 2g: BLAST + powerneedle wGRR per genome vs reference K-locus; writes `locus/grr.tsv` |
| `gwas_sgnh_selection.py` | `select_sgnh_predictors()` | Step 3: best SGNH predictor per K-locus; writes `gwas_sgnh_hits.tsv` |

## Run

```
cd scripts/processing/gwas-proc
conda run -n jkoszucki python main.py
```
