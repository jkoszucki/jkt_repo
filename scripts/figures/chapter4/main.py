from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from branch_core_heatmap import plot_branch_core_heatmap
from composition_groups import plot_composition_groups
from config import Config
from cytoscape_export import export_cytoscape
from jp_heatmap import plot_jp_heatmap
from jp_scatter import plot_jp_scatter
from ktypes_draw import KTypeStructureDrawer
from ktypes_plots import KTypePlotAPI
from figure2_table import plot_figure2_table
from mod_freq import plot_modification_frequency
from mw_hist import plot_mw_hist


def main() -> None:
    cfg = Config()
    chapter3_output_dir = cfg.output_dir / "chapter3"
    chapter2_dir = Path(__file__).resolve().parent

    processed_csv = chapter3_output_dir / "ktypes.csv"
    similarity_csv = chapter3_output_dir / "ktypes_sim.csv"
    modifications_csv = chapter3_output_dir / "ktypes_modifications.csv"

    plots_dir = chapter2_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    ### DRAW K-TYPE STRUCTURES (checkpoint: skip already generated)
    drawer = KTypeStructureDrawer(
        output_dir=chapter2_dir,
        processed_csv=processed_csv,
        modifications_csv=modifications_csv,
        plots_subdir="plots/ktypes_structures",
    )
    ktypes_df = pd.read_csv(processed_csv)
    skipped = []
    for _, row in ktypes_df.iterrows():
        kname = str(row.get("structure_id", "") or row.get("K_type", "")).strip()
        stem = (kname if kname else f"row_{_}").replace("/", "_").replace(" ", "_")
        if (drawer.png_dir / f"{stem}.png").exists() and (drawer.svg_dir / f"{stem}.svg").exists():
            skipped.append(kname)
        else:
            drawer.draw_single(row, stem_name=stem)
    if skipped:
        print(f"Skipped {len(skipped)} already-generated structures: {', '.join(skipped)}")

    ### PLOTS
    plot_api = KTypePlotAPI(
        ktypes_csv=processed_csv,
        similarity_csv=similarity_csv,
        modifications_csv=modifications_csv,
        style=cfg.style,
    )

    plot_api.plot_cumulative_structures(
        output_path=plots_dir / "year.png", show=False
    )

    plot_api.plot_similarity_histogram(
        output_path=plots_dir / "score_hist.png", column="path_jaccard_total", bins=40, show=False
    )

    plot_api.plot_modification_and_monosaccharide_distribution(
        output_path=plots_dir / "mono_prevalence.png",
        show=False,
    )

    ### QUICK PLOTS
    plot_branch_core_heatmap(
        similarity_csv=similarity_csv,
        output_path=plots_dir / "branch_core_heatmap.png",
        style=cfg.style,
    )

    plot_composition_groups(
        input_csv=processed_csv,
        outdir=plots_dir / "composition",
    )

    plot_mw_hist(
        input_csv=processed_csv,
        output_path=plots_dir / "mw_hist.png",
        style=cfg.style,
    )

    plot_modification_frequency(
        modifications_csv=modifications_csv,
        output_path=plots_dir / "modification_freq.png",
        style=cfg.style,
    )

    plot_jp_scatter(
        similarity_csv=similarity_csv,
        output_path=plots_dir / "jp_scatter.png",
        style=cfg.style,
    )

    plot_jp_heatmap(
        similarity_csv=similarity_csv,
        output_path=plots_dir / "jp_heatmap.png",
        style=cfg.style,
    )

    export_cytoscape(
        similarity_csv=similarity_csv,
        ktypes_csv=processed_csv,
        modifications_csv=modifications_csv,
        output_dir=plots_dir / "cytoscape",
        style=cfg.style,
    )

    plot_figure2_table(
        input_xlsx=cfg.input_dir / "cps.xlsx",
        output_path=plots_dir / "figure2_table.png",
        style=cfg.style,
    )


if __name__ == "__main__":
    main()
