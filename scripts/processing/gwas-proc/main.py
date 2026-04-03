"""
processing/gwas-proc — shared data preparation for chapters 2 and 3.

Pipeline:
  1. filter_and_classify()           → gwas_hits.tsv
  2. per-PC file export
       2a. export_per_locus()            → {ecod_folder}/{locus}/{cl}/{PC}/protein/pc.fasta, structure.cif
       2b. export_per_protein()          → protein/sequence.fasta
       2c. prepare_af3_batches()         → manual-upload/gwas_batch_NNN.json  [optional: RUN_AF3_JSON_PREP]
       2d. prepare_prophage_db()         → manual-outputs/prophage/prophage_proteins.faa, prophage_db.*
       2e. blast_repr_vs_prophage()      → protein/against-prophages/raw_blast.tsv
       2f. export_klocus_genbank_dataset() → locus/reference.gb, match/, locus_only/, pc_only/
       2g. export_grr_vs_reference()     → locus/grr.tsv  [optional: RUN_GRR_VS_REFERENCE]
  3. select_sgnh_predictors()        → gwas_sgnh_hits.tsv

Output root: cfg.output_dir / "gwas-data"
Feeds figures/chapter2/ (no-ecod, ssrbh-ecod entries)
and  figures/chapter3/ (sgnh-ecod entries).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from gwas_filter_classify import filter_and_classify
from per_locus_export import (
    export_per_locus,
    export_per_protein,
    blast_repr_vs_prophage,
    export_klocus_genbank_dataset,
)
from per_locus_grr import export_grr_vs_reference
from gwas_sgnh_selection import select_sgnh_predictors
from prophage_db import prepare_prophage_db
from af3_json_prep import prepare_af3_batches

cfg = Config()

# ---------------------------------------------------------------------------
# Run flags — toggle steps without modifying config.yml
# When ON: step runs with per-PC checkpointing — skips any PC whose output
# file already exists (raw_blast.tsv / grr.tsv). Delete files to recompute.
# When OFF: step is skipped entirely.
# ---------------------------------------------------------------------------
RUN_AF3_JSON_PREP    = True    # True → write gwas_batch_NNN.json; skips PCs where structure.cif already exists
RUN_PROPHAGE_BLAST   = True    # True → blastp per PC; skips if raw_blast.tsv exists
RUN_GRR_VS_REFERENCE = True    # True → BLAST + powerneedle per PC; skips if grr.tsv exists

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
table_path     = cfg.input_dir / "gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv"
functions_path = cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/3_FUNCTIONS/clusters_functions_best_all.tsv"
mmseqs_dir     = cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/2_MMSEQS"
gwas_root      = cfg.output_dir / "gwas-data"

af3_dir      = cfg.output_dir / "manual-outputs" / "alphafold3"
blast_db_dir        = cfg.output_dir / "manual-outputs" / "prophage"

# ---------------------------------------------------------------------------
# Step 1 — filter, classify, write gwas_hits.tsv
# ---------------------------------------------------------------------------
print("Step 1: Filtering and classifying GWAS hits …")
filter_and_classify(table_path, functions_path, gwas_root)

gwas_hits_tsv = gwas_root / "gwas_hits.tsv"

# ---------------------------------------------------------------------------
# Step 2 — per-PC file export
# ---------------------------------------------------------------------------
print("\nStep 2a: Exporting per-PC alignment FASTAs and AF3 structures …")
export_per_locus(
    gwas_hits_tsv = gwas_hits_tsv,
    mmseqs_dir    = mmseqs_dir,
    af3_dir       = af3_dir,
    gwas_root     = gwas_root,
)

print("\nStep 2b: Writing representative sequence FASTAs …")
export_per_protein(
    gwas_hits_tsv = gwas_hits_tsv,
    mmseqs_dir    = mmseqs_dir,
    gwas_root     = gwas_root,
)

if RUN_AF3_JSON_PREP:
    print("\nStep 2c: Preparing AF3 JSON upload files …")
    sequence_paths = sorted(gwas_root.rglob("protein/sequence.fasta"))
    pending = [p for p in sequence_paths if not (p.parent / "structure.cif").exists()]
    print(f"  {len(sequence_paths)} sequences found; {len(pending)} pending AF3 upload (no structure.cif yet)")
    prepare_af3_batches(
        sequence_paths = pending,
        output_dir     = cfg.output_dir,
    )
else:
    print("\n[skip] Step 2c: AF3 JSON prep (RUN_AF3_JSON_PREP = False)")

print("\nStep 2d: Building prophage BLAST DB …")
prepare_prophage_db(
    prophage_faa_pattern = str(cfg.input_dir / "gwas/2_PROPHAGES/4_FASTA_CDS_AA/*.faa"),
    db_dir               = blast_db_dir,
)

print("\nStep 2e: BLASTP vs prophage proteins …")
blast_repr_vs_prophage(
    gwas_hits_tsv        = gwas_hits_tsv,
    mmseqs_dir           = mmseqs_dir,
    prophage_faa_pattern = str(cfg.input_dir / "gwas/2_PROPHAGES/4_FASTA_CDS_AA/*.faa"),
    gwas_root            = gwas_root,
    blast_db_dir         = blast_db_dir,
    run_blast            = RUN_PROPHAGE_BLAST,
)

print("\nStep 2f: Exporting K-locus GenBank dataset …")
export_klocus_genbank_dataset(
    gwas_hits_tsv         = gwas_hits_tsv,
    bacteria_metadata_tsv = cfg.input_dir / "gwas/1_BACTERIA/bacteria_metadata.tsv",
    k_loci_dir            = cfg.input_dir / "gwas/4_K_LOCI",
    ref_dir               = cfg.input_dir / "cps/cps_loci_ref/K_locus_primary_referenece_MGG",
    mmseqs_dir            = mmseqs_dir,
    gwas_root             = gwas_root,
)

if RUN_GRR_VS_REFERENCE:
    print("\nStep 2g: Computing wGRR vs reference K-locus …")
    export_grr_vs_reference(
        gwas_hits_tsv = gwas_hits_tsv,
        gwas_root     = gwas_root,
    )
else:
    print("\n[skip] Step 2g: wGRR vs reference (RUN_GRR_VS_REFERENCE = False)")

# ---------------------------------------------------------------------------
# Step 3 — Best SGNH predictor per K-locus
# ---------------------------------------------------------------------------
print("\nStep 3: Selecting best SGNH predictors per K-locus …")
select_sgnh_predictors(gwas_hits_tsv, mmseqs_dir, gwas_root)

print("\nDone.")
