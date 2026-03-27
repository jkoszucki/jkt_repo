from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from sgnh_recall_plot import plot_sgnh_recall
from sgnh_network import compute_sgnh_network
from phandango_export import export_phandango
from phandango_legend import plot_phandango_legend
from phandango_render import render_all_phandango
from structure_render import render_all_structures, render_af3_no_ecod, render_af3_no_ecod_plddt
from functions_heatmap import plot_functions_heatmap
from alignment_profiles import plot_alignment_profiles
from grr_vs_reference_plot import plot_grr_vs_reference_all, plot_grr_vs_reference_global
from gwas_scatter_plot import plot_gwas_scatter, plot_gwas_phrog_barplot

# ---------------------------------------------------------------------------
# Domain colour configuration for AF3 structure rendering
#
# Each entry: (residue_selection, color)
#   residue_selection — PyMOL resi syntax, e.g. "1-200" or "1-50+100-150"
#   color             — hex string, or None → uses style.sgnh_domain_color
#
# Unlisted residues keep the default gray (#cfcfce).
# Leave as an empty list to render the whole structure in gray.
# ---------------------------------------------------------------------------
DOMAIN_COLORS: list[tuple[str, str | None]] = [
    # ("1-350", None),   # SGNH hydrolase domain → sgnh_domain_color from config
]

# ---------------------------------------------------------------------------
# Run flags — toggle steps without modifying config.yml
# When ON: always executes and overwrites existing files.
# ---------------------------------------------------------------------------
NO_ECOD_PCS = [
    "PC0232", "PC0294", "PC0370", "PC0445", "PC0467", "PC0478",
    "PC0585", "PC0597", "PC0666", "PC0712", "PC1259", "PC1260",
]
AF3_FOLDS_DIR = "chapter2/raw_af3/folds_2026_03_26_10_27"

PLOT_ALIGNMENT_PROFILES = True
RENDER_PHANDANGO        = True
PLOT_GRR_VS_REFERENCE   = True
PLOT_GWAS_SCATTER       = True


def main() -> None:
    cfg = Config()
    chapter2_output_dir = cfg.output_dir / "chapter2"
    plots_dir = Path(__file__).resolve().parent / "plots"
    plots_dir.mkdir(exist_ok=True)

    best_predictors_csv = chapter2_output_dir / "best_predictors.csv"

    experimental_fastas = [
        cfg.input_dir / "gwas/sgnh-active.fasta",
    ]

    plot_sgnh_recall(best_predictors_csv, plots_dir, style=cfg.style)

    if PLOT_ALIGNMENT_PROFILES:
        print("Plotting alignment gap-frequency profiles …")
        plot_alignment_profiles(
            per_locus_dir       = chapter2_output_dir / "per-locus",
            best_predictors_csv = best_predictors_csv,
            plots_dir           = plots_dir,
            style               = cfg.style,
        )
    else:
        print("[skip] alignment profiles (PLOT_ALIGNMENT_PROFILES = False)")

    plot_functions_heatmap(
        functions_tsv       = chapter2_output_dir / "best_predictors_functions.tsv",
        best_predictors_csv = best_predictors_csv,
        plots_dir           = plots_dir,
        style               = cfg.style,
    )
    compute_sgnh_network(
        best_predictors_csv = best_predictors_csv,
        experimental_fastas = experimental_fastas,
        output_dir          = plots_dir,
        tmp_dir             = plots_dir / "sgnh_network_blast_tmp",
    )

    plot_phandango_legend(plots_dir, style=cfg.style)

    print("Exporting phandango files …")
    export_phandango(
        best_predictors_csv = best_predictors_csv,
        phandango_src_root  = cfg.input_dir / "gwas/3_GWAS/4_ANALYZE/1_PER_LOCUS",
        mmseqs_dir          = cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/2_MMSEQS",
        plots_dir           = plots_dir,
        style               = cfg.style,
    )

    if RENDER_PHANDANGO:
        print("Rendering phandango PNGs …")
        render_all_phandango(
            best_predictors_csv = best_predictors_csv,
            plots_dir           = plots_dir,
            meta_tsv            = cfg.input_dir / "gwas/1_BACTERIA/bacteria_metadata.tsv",
            style               = cfg.style,
        )
    else:
        print("[skip] phandango PNGs (RENDER_PHANDANGO = False)")

    if PLOT_GRR_VS_REFERENCE:
        print("Plotting GRR vs reference …")
        plot_grr_vs_reference_all(
            best_predictors_csv = best_predictors_csv,
            per_locus_dir       = chapter2_output_dir / "per-locus",
            plots_dir           = plots_dir,
            style               = cfg.style,
        )
        plot_grr_vs_reference_global(
            best_predictors_csv = best_predictors_csv,
            per_locus_dir       = chapter2_output_dir / "per-locus",
            plots_dir           = plots_dir,
            style               = cfg.style,
        )
    else:
        print("[skip] GRR vs reference plots (PLOT_GRR_VS_REFERENCE = False)")

    if PLOT_GWAS_SCATTER:
        print("Plotting GWAS scatter (F1 vs MCC) …")
        plot_gwas_scatter(
            pyseer_hits_tsv = cfg.input_dir / "gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv",
            hhsearch_dir    = cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/3_FUNCTIONS/HHSEARCH",
            plots_dir       = plots_dir,
            style           = cfg.style,
        )
        plot_gwas_phrog_barplot(
            pyseer_hits_tsv = cfg.input_dir / "gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv",
            hhsearch_dir    = cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/3_FUNCTIONS/HHSEARCH",
            plots_dir       = plots_dir,
            style           = cfg.style,
        )
    else:
        print("[skip] GWAS scatter (PLOT_GWAS_SCATTER = False)")

    print("Rendering structure PNGs …")
    render_all_structures(
        per_locus_dir = chapter2_output_dir / "per-locus",
        plots_dir     = plots_dir,
        domain_colors = DOMAIN_COLORS,
        style         = cfg.style,
    )

    print("Rendering AF3 no-ECOD structure PNGs …")
    render_af3_no_ecod(
        af3_folds_dir = cfg.output_dir / AF3_FOLDS_DIR,
        pcs           = NO_ECOD_PCS,
        plots_dir     = plots_dir,
        style         = cfg.style,
    )

    print("Rendering AF3 no-ECOD pLDDT-coloured PNGs …")
    render_af3_no_ecod_plddt(
        af3_folds_dir = cfg.output_dir / AF3_FOLDS_DIR,
        pcs           = NO_ECOD_PCS,
        plots_dir     = plots_dir,
        style         = cfg.style,
    )


if __name__ == "__main__":
    main()
