# processing/acetyl-proc — acetyltransferase detection and annotation

Chapter 2 specific data preparation. Detects acetyltransferases at two levels (prophage PC80 clusters and GWAS predictors), exports per-PC annotation, and runs structural comparisons for no-ecod GWAS predictors. Feeds `figures/chapter2/`.

## Pipeline

```
raw_hhsuite.tsv  +  pc2proteins.tsv
        │
        ▼ Step 1: SSLBH/acetyltransferase detection at PC80 level (three evidence sources)
        │   ECOD:   db == 'ECOD',   name contains 'Single-stranded left-handed beta-helix', prob ≥ 0.7
        │   PHROGs: db == 'PHROGS', annot in {acetyltransferase, O-acetyltransferase,
        │           O-antigen acetylase, acyl-CoA N-acyltransferase,
        │           antimicrobial peptide resistance and lipid A acylation protein PagP,
        │           FabG-like 3-oxoacyl-(acyl-carrier-protein) reductase}, prob ≥ 0.9
        │   PFAM:   db == 'PFAM',   name contains 'PF01757', prob ≥ 0.9
        │   Map detected PC80s → protein IDs via pc2proteins
        │   202 PC80s detected → 2,791 SSLBH protein IDs
        │
        ▼ Step 2: flag pyseer GWAS hits with SSLBH evidence (HHsearch)
        │   Load pyseer_hits_all.tsv; filter pvalue_corr ≤ 0.05
        │   For each clustering level: scan alignment FASTAs → find GWAS PCs containing SSLBH proteins
        │   Add columns: sslbh_detected (bool), sslbh_pc80, sslbh_detected_by
        │   → acetyltransferase/acetyl-gwas/pyseer_hits_sslbh_pvalcor005.tsv
        │
        ▼ Step 3: PC80 annotation export
        │   → acetyltransferase/acetyl-annot-pc80/acetyl-hhsearch-pc80.tsv  (raw_hhsuite rows for detected PC80s)
        │   → acetyltransferase/acetyl-annot-pc80/acetyl-pc80.tsv      (pc80 → proteinID → prophageID → genomeID → sequence)
        │   → acetyltransferase/acetyl-annot-pc80/{pc80}/
        │         sequence.fasta, sequence.cif (if available),
        │         against-prophages/raw_blast.tsv
        │
        ▼ Step 4: TM-align screen for no-ecod GWAS predictors
        │   Query: AF3 monomer for each no-ecod GWAS PC (uploaded manually → structure.cif)
        │   Reference: PROTEIN01_MOD_AC_K1, PROTEIN02_MOD_AC_K2, PROTEIN03_MOD_AC_K57
        │   Invoked via subprocess: conda run -n tmtools python lib/tmalign_screen.py
        │   → acetyltransferase/acetyl-gwas/no-ecod-reported-topology-tm-align.tsv
        │
        └─▶ Step 5: merge TM-align results into pyseer table
                Read no-ecod-reported-topology-tm-align.tsv; deduplicate by (clustering_level, PC, reference)
                Pivot wide: one row per (clustering_level, PC) with columns per reference
                For PCs with any tm_score_query > 0: sslbh_detected = True, sslbh_detected_by += "TM-align"
                Adds {reference}_tm_score_query/tm_score_ref/rmsd columns to pyseer table
                → rewrites acetyltransferase/acetyl-gwas/pyseer_hits_sslbh_pvalcor005.tsv
```

## Run flags

Defined at the top of `main.py`:

| Flag | Default | Controls |
|------|---------|----------|
| `RUN_PROPHAGE_BLAST` | `True` | BLASTP per PC80 representative vs prophage DB; skips if `raw_blast.tsv` exists |
| `RUN_TMALIGN` | `True` | TM-align screen for no-ecod GWAS predictors (Step 4); skips individual rows already in `no-ecod-reported-topology-tm-align.tsv`; Step 5 merge always runs if file exists |

## Input data

| Path | Description |
|------|-------------|
| `input_dir/gwas/2_PROPHAGES/raw_hhsuite.tsv` | HHsearch results at PC80 level; key columns: `query` (PC80), `db`, `name`, `annot`, `prob` |
| `input_dir/gwas/2_PROPHAGES/pcs2proteins.tsv` | PC80 → protein ID mapping; columns: `PC`, `proteinID`, `repr` |
| `input_dir/gwas/2_PROPHAGES/prophages_metadata.tsv` | columns `prophageID`, `genomeID` — maps protein IDs to genomes |
| `input_dir/gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv` | All pyseer GWAS hits; `version` column = clustering level |
| `input_dir/gwas/3_GWAS/1_INTERMEDIATE/2_MMSEQS/{clustering_level}/alignments/` | Per-PC alignment FASTAs for SSLBH protein → GWAS PC mapping |
| `input_dir/gwas/2_PROPHAGES/4_FASTA_CDS_AA/*.faa` | Prophage protein FASTAs (shared BLAST DB; protein sequences for PC80 export) |
| `output_dir/gwas-data/gwas_hits.tsv` | GWAS hits passing quality filters (used by TM-align step to identify no-ecod predictors) |
| `output_dir/gwas-data/no-ecod-reported-topology/` | AF3 structures for no-ecod GWAS predictors (TM-align query) |
| `output_dir/manual-outputs/alphafold3/experimental/` | AF3 monomers for reference acetyltransferases (TM-align references) |

## Output data

Written to `output_dir/acetyltransferase/`.

```
acetyltransferase/
    acetyl-gwas/
        pyseer_hits_sslbh_pvalcor005.tsv   — all pyseer rows with pvalue_corr ≤ 0.05;
                                              sslbh_detected (bool), sslbh_pc80, sslbh_detected_by
                                              (from HHsearch PC80 detection and/or TM-align);
                                              per-reference TM-align columns appended by Step 5:
                                              {reference}_tm_score_query, {reference}_tm_score_ref,
                                              {reference}_rmsd  (NaN for non-no-ecod rows)

        no-ecod-reported-topology-tm-align.tsv — raw TM-align results; one row per
                                              (locus, clustering_level, pc, reference);
                                              columns: locus, clustering_level, pc, reference,
                                                       tm_score_query, tm_score_ref, rmsd

    acetyl-annot-pc80/
        acetyl-hhsearch-pc80.tsv           — raw_hhsuite.tsv rows for all 202 detected PC80s
        acetyl-pc80.tsv                    — one row per protein:
                                              pc80, detected_by, proteinID, prophageID,
                                              genomeID, sequence
        {pc80}/
            sequence.fasta                 — representative protein (first in pc2proteins)
            sequence.cif                   — AF3 monomer if found in manual-outputs/alphafold3/
            against-prophages/
                raw_blast.tsv              — BLASTP vs prophage DB; columns:
                                              qseqid, sseqid, pident, length, mismatch,
                                              gapopen, qstart, qend, sstart, send,
                                              evalue, bitscore, qcovhsp, qlen, slen, sseq
```

## Detection criteria

### Step 1 — PC80 detection

All three sources are queried from `raw_hhsuite.tsv`. A PC80 is included if it passes **any** filter. The `detected_by` value records all sources that fired (comma-separated, e.g. `ECOD`, `PHROGs,PFAM`).

| Source | Column | Filter |
|--------|--------|--------|
| ECOD | `name` | `db == 'ECOD'` AND `name` contains `Single-stranded left-handed beta-helix` AND `prob ≥ 0.7` |
| PHROGs | `annot` | `db == 'PHROGS'` AND `annot` (case-insensitive) in `{acetyltransferase, o-acetyltransferase, o-antigen acetylase, acyl-CoA N-acyltransferase, antimicrobial peptide resistance and lipid A acylation protein PagP, FabG-like 3-oxoacyl-(acyl-carrier-protein) reductase}` AND `prob ≥ 0.9` |
| PFAM | `name` | `db == 'PFAM'` AND `name` contains `PF01757` AND `prob ≥ 0.9` |

Note: PHROGs functional annotation is in the `annot` column, not `name`. The `name` column for PHROGs contains the profile identifier (e.g. `phrog_1249 ## NC_028770_p48`).

Protein IDs for detected PC80s are pulled from `pc2proteins.tsv` (`PC` and `proteinID` columns).

### Step 2 — GWAS flagging

```
proteinID  →  strip _PROTEIN_\d+ suffix  →  prophageID
prophageID  →  look up in prophages_metadata.tsv  →  genomeID
```

For each clustering level, alignment FASTAs (`{clustering_level}/alignments/{PC}.fasta`) are scanned by header only for efficiency. A GWAS PC is SSLBH-positive if at least one of its member proteins belongs to a detected PC80.

Added columns in `pyseer_hits_sslbh_pvalcor005.tsv`:

| Column | Description |
|--------|-------------|
| `sslbh_detected` | `True` if the GWAS PC contains ≥1 SSLBH protein |
| `sslbh_pc80` | comma-separated PC80 IDs responsible for the hit |
| `sslbh_detected_by` | comma-separated sources across matched PC80s |

### Step 4 — TM-align screen

Applies only to GWAS predictors with `ecod_type == "no-ecod"` from `gwas_hits.tsv`, whose AF3 monomer (`structure.cif`) has been placed at:
```
output_dir/gwas-data/no-ecod-reported-topology/{locus}/{clustering_level}/{pc}/protein/structure.cif
```

`main.py` invokes `tmalign_screen.py` via subprocess in the `tmtools` conda environment:

```python
subprocess.run(
    ["conda", "run", "-n", "tmtools", "python",
     str(Path(__file__).parent / "lib" / "tmalign_screen.py"),
     "--gwas-hits",        str(gwas_hits_tsv),
     "--gwas-root",        str(gwas_root),
     "--experimental-dir", str(experimental_dir),
     "--output",           str(acetyl_dir / "acetyl-gwas" / "no-ecod-reported-topology-tm-align.tsv")],
    check=True,
)
```

Reference proteins: `PROTEIN01_MOD_AC_K1`, `PROTEIN02_MOD_AC_K2`, `PROTEIN03_MOD_AC_K57` at `output_dir/manual-outputs/alphafold3/enzymes/{lowercase_protein_id}/fold_{lowercase_protein_id}_model_0.cif`. Both TM-scores (query-normalised and reference-normalised) and RMSD are written. No threshold is applied — the full results table is the output.

If reference structures are absent, the step is skipped with a warning and `main.py` exits normally.

### Step 5 — TM-align merge into pyseer table

Reads `acetyl-gwas/no-ecod-reported-topology-tm-align.tsv`. Deduplicates by `(clustering_level, PC, reference)` (same structure → same score), then pivots wide to one row per `(clustering_level, PC)`. Merged into `pyseer_hits_sslbh_pvalcor005.tsv` via left join on `(clustering_level, PC)`.

For each PC with `tm_score_query > 0` for any reference:
- `sslbh_detected` is set to `True`
- `"TM-align"` is added to `sslbh_detected_by` (comma-separated, deduplicated)
- `sslbh_pc80` is left empty (not applicable — structural detection, not PC80 HHsearch)

If `no-ecod-reported-topology-tm-align.tsv` does not exist (references not yet placed), Step 5 is skipped with a warning; the pyseer table is unchanged.

## Clustering levels

| Level | Identity | Coverage |
|-------|----------|----------|
| `PCI00C50` | 0% | 50% |
| `PCI50C50` | 50% | 50% |
| `PCI80C50` | 80% | 50% |
| `PCI00C80` | 0% | 80% |
| `PCI50C80` | 50% | 80% |
| `PCI80C80` | 80% | 80% |

PC80 refers to the `PCI80C50` level used by `raw_hhsuite.tsv` and `pcs2proteins.tsv`.

## Library modules (`lib/`)

| Module | Function | Role |
|--------|----------|------|
| `sslbh_detect.py` | `detect_sslbh_pc80()` | Step 1: query raw_hhsuite for ECOD/PHROGs/PFAM hits; returns DataFrame with pc80, detected_by, protein_ids list |
| `cluster_map.py` | `flag_pyseer_with_sslbh()` | Step 2: scan FASTAs per clustering level, flag SSLBH-positive rows in pyseer table, write output TSV |
| `pc80_export.py` | `export_pc80_annotation()` | Step 3: write acetyl-hhsearch-pc80.tsv, acetyl-pc80.tsv, per-PC sequence/BLAST files |
| `tmalign_screen.py` | CLI script | Step 4: TM-align no-ecod GWAS structures vs reference enzymes; invoked via subprocess in `tmtools` env |
| `merge_tmalign.py` | `merge_tmalign_into_pyseer()` | Step 5: pivot TM-align results wide, merge into pyseer table, update sslbh_detected/sslbh_detected_by |
| `helpers/prophage_db.py` | `prepare_prophage_db()` | Build shared BLAST DB and sequence index (shared with gwas-proc) |

## Run

```
cd scripts/processing/acetyl-proc
conda run -n jkoszucki python main.py
```

Step 4 (TM-align) is invoked automatically via subprocess in the `tmtools` environment. Re-run after placing reference structures to complete the TM-align screen.
