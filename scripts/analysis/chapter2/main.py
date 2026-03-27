import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from gwas_predictor_selection import select_best_predictors
from per_locus_export import export_per_locus, export_klocus_genbank_dataset, blast_repr_vs_prophage, export_per_protein
from per_locus_grr import export_grr_vs_reference

cfg = Config()

# ---------------------------------------------------------------------------
# Run flags — toggle steps without modifying config.yml
# ---------------------------------------------------------------------------
RUN_PROPHAGE_BLAST      = False
RUN_GRR_VS_REFERENCE    = True   # True → runs BLAST + powerneedle per locus, overwrites grr_vs_reference.tsv

table_path   = cfg.input_dir / "gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv"
mmseqs_dir   = cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/2_MMSEQS"
hhsearch_dir = cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/3_FUNCTIONS/HHSEARCH"
output_dir   = cfg.output_dir / "chapter2"

select_best_predictors(table_path, mmseqs_dir, hhsearch_dir, output_dir)

print("\nExporting per-locus data ...")
export_per_locus(
    best_predictors_csv = output_dir / "best_predictors.csv",
    mmseqs_dir          = mmseqs_dir,
    raw_af3_dir         = output_dir / "raw_af3",
    experimental_fastas = [cfg.input_dir / "gwas/sgnh-active.fasta"],
    output_dir          = output_dir / "per-locus",
)

experimental_fastas = [cfg.input_dir / "gwas/sgnh-active.fasta"]

print("\nBLASTing representative sequences vs prophage proteins ...")
blast_repr_vs_prophage(
    best_predictors_csv  = output_dir / "best_predictors.csv",
    prophage_faa_pattern = str(cfg.input_dir / "gwas/2_PROPHAGES/4_FASTA_CDS_AA/*.faa"),
    per_locus_dir        = output_dir / "per-locus",
    experimental_fastas  = experimental_fastas,
    blast_db_dir         = Path("/Users/januszkoszucki/Projects/thesis-claude/data/playground/prophage_blast_db"),
    run_blast            = RUN_PROPHAGE_BLAST,
)

print("\nExporting per-protein data ...")
export_per_protein(
    best_predictors_csv = output_dir / "best_predictors.csv",
    per_locus_dir       = output_dir / "per-locus",
    experimental_fastas = experimental_fastas,
)

print("\nExporting per-locus K-locus GenBank dataset ...")
export_klocus_genbank_dataset(
    best_predictors_csv   = output_dir / "best_predictors.csv",
    bacteria_metadata_tsv = cfg.input_dir / "gwas/1_BACTERIA/bacteria_metadata.tsv",
    k_loci_dir            = cfg.input_dir / "gwas/4_K_LOCI",
    ref_dir               = cfg.input_dir / "cps/cps_loci_ref/K_locus_primary_referenece_MGG",
    mmseqs_dir            = mmseqs_dir,
    output_dir            = output_dir / "per-locus",
)

if RUN_GRR_VS_REFERENCE:
    print("\nComputing wGRR vs reference per locus ...")
    export_grr_vs_reference(
        best_predictors_csv = output_dir / "best_predictors.csv",
        per_locus_dir       = output_dir / "per-locus",
    )
else:
    print("\n[skip] GRR vs reference (RUN_GRR_VS_REFERENCE = False)")
