from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from gwas_scatter_plot import plot_gwas_scatter, plot_gwas_phrog_barplot
from no_ecod_wordcloud import plot_no_ecod_wordcloud
from no_ecod_network import compute_no_ecod_network
from sgnh_recall_plot import plot_sgnh_recall
from sgnh_network import compute_sgnh_network
from phandango_export import export_phandango
from phandango_legend import plot_phandango_legend
from phandango_render import render_all_phandango
from structure_renderer import StructureRenderer
from functions_heatmap import plot_functions_heatmap
from alignment_profiles import plot_alignment_profiles
from grr_vs_reference_plot import plot_grr_vs_reference_all, plot_grr_vs_reference_global

# ---------------------------------------------------------------------------
# Run flags — toggle steps without modifying config.yml
# When ON: always executes and overwrites existing files.
# ---------------------------------------------------------------------------
RENDER_MODE = "fast"    # "fast" → no ray tracing; "accurate" → ray tracing (slow)

# Chapter 2 — Figure 1
PLOT_FIGURE1_PANELA     = True    # 6-panel F1 vs MCC scatter per clustering level
PLOT_FIGURE1_PANELB     = False   # word cloud — requires FoldSeek TSV in output_dir/raw/
PLOT_FIGURE1_PANELC     = False   # no-ecod similarity network — requires clustered FASTA

# Chapter 2 — exploratory / legacy
PLOT_GWAS_PHROG_BAR     = False   # PHROG function bar plot (PCI50C50 only)

# Chapter 3 content (temporarily lives here; will move to figures/chapter3/)
PLOT_SGNH_RECALL        = True
COMPUTE_SGNH_NETWORK    = True
PLOT_FUNCTIONS_HEATMAP  = True
EXPORT_PHANDANGO        = True
RENDER_PHANDANGO        = True
PLOT_ALIGNMENT_PROFILES = True
PLOT_GRR_VS_REFERENCE   = True

# Structure rendering (requires pymol conda env)
RENDER_STRUCTURES       = False


def main() -> None:
    cfg = Config()
    chapter2_output_dir = cfg.output_dir / "chapter2"
    plots_dir = Path(__file__).resolve().parent / "plots"
    plots_dir.mkdir(exist_ok=True)

    best_predictors_csv = chapter2_output_dir / "best_predictors.csv"
    experimental_fastas = [cfg.input_dir / "gwas/sgnh-active.fasta"]

    # -----------------------------------------------------------------------
    # Chapter 2 — Figure 1
    # -----------------------------------------------------------------------

    if PLOT_FIGURE1_PANELA:
        print("Figure 1 Panel A: 6-panel GWAS scatter …")
        plot_gwas_scatter(
            pyseer_hits_tsv = cfg.input_dir / "gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv",
            hhsearch_dir    = cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/3_FUNCTIONS/HHSEARCH",
            plots_dir       = plots_dir,
            style           = cfg.style,
        )
    else:
        print("[skip] Figure 1 Panel A (PLOT_FIGURE1_PANELA = False)")

    if PLOT_FIGURE1_PANELB:
        print("Figure 1 Panel B: no-ecod word cloud …")
        plot_no_ecod_wordcloud(
            foldseek_tsv = chapter2_output_dir / "raw" / "no_ecod_foldseek.tsv",
            plots_dir    = plots_dir,
            style        = cfg.style,
        )
    else:
        print("[skip] Figure 1 Panel B (PLOT_FIGURE1_PANELB = False) "
              "— place FoldSeek TSV at output_dir/chapter2/raw/no_ecod_foldseek.tsv")

    if PLOT_FIGURE1_PANELC:
        print("Figure 1 Panel C: no-ecod similarity network …")
        compute_no_ecod_network(
            sequences_fasta   = chapter2_output_dir / "raw" / "no_ecod_representatives.fasta",
            plots_dir         = plots_dir,
            tmp_dir           = plots_dir / "no_ecod_network_blast_tmp",
            node_metadata_csv = chapter2_output_dir / "raw" / "no_ecod_metadata.csv",
        )
    else:
        print("[skip] Figure 1 Panel C (PLOT_FIGURE1_PANELC = False) "
              "— place clustered FASTA at output_dir/chapter2/raw/no_ecod_representatives.fasta")

    if PLOT_GWAS_PHROG_BAR:
        print("GWAS PHROG functions bar plot …")
        plot_gwas_phrog_barplot(
            pyseer_hits_tsv = cfg.input_dir / "gwas/3_GWAS/3_PROCESSING/pyseer_hits_all.tsv",
            hhsearch_dir    = cfg.input_dir / "gwas/3_GWAS/1_INTERMEDIATE/3_FUNCTIONS/HHSEARCH",
            plots_dir       = plots_dir,
            style           = cfg.style,
        )

    # -----------------------------------------------------------------------
    # Chapter 3 content — temporarily here; will move to figures/chapter3/
    # -----------------------------------------------------------------------

    if PLOT_SGNH_RECALL:
        plot_sgnh_recall(best_predictors_csv, plots_dir, style=cfg.style)

    if PLOT_FUNCTIONS_HEATMAP:
        plot_functions_heatmap(
            functions_tsv       = chapter2_output_dir / "best_predictors_functions.tsv",
            best_predictors_csv = best_predictors_csv,
            plots_dir           = plots_dir,
            style               = cfg.style,
        )

    if COMPUTE_SGNH_NETWORK:
        compute_sgnh_network(
            best_predictors_csv = best_predictors_csv,
            experimental_fastas = experimental_fastas,
            output_dir          = plots_dir,
            tmp_dir             = plots_dir / "sgnh_network_blast_tmp",
        )

    plot_phandango_legend(plots_dir, style=cfg.style)

    if EXPORT_PHANDANGO:
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

    # -----------------------------------------------------------------------
    # Structure rendering — requires pymol conda env
    # -----------------------------------------------------------------------

    if RENDER_STRUCTURES:
        renderer = StructureRenderer(render_mode=RENDER_MODE, style=cfg.style)

        print("Rendering GWAS structure PNGs …")
        gwas_root = cfg.output_dir / "gwas-data"
        cif_paths = sorted(gwas_root.rglob("*/protein/structure.cif"))
        renderer.render_next_to(cif_paths)

        print("Rendering acetyltransferase structure PNGs …")
        acetyl_cifs = sorted((cfg.output_dir / "acetyltransferase").rglob("structure.cif"))
        renderer.render_next_to(acetyl_cifs)


if __name__ == "__main__":
    main()
