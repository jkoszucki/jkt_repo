from __future__ import annotations

from pathlib import Path

from ktypes_draw import KTypeStructureDrawer
from ktypes_modifications import KTypeModificationsAPI
from ktypes_process import KTypeTablesAPI
from ktypes_similarity import KTypeSimilarityAPI
from plot_ktypes_api import KTypePlotAPI

ANALYSIS_DIR = (Path(__file__).resolve().parent / "../analysis").resolve()
INPUT_XLSX = ANALYSIS_DIR / "ktypes.xlsx"
PROCESSED_CSV = ANALYSIS_DIR / "ktypes.csv"
SIMILARITY_CSV = ANALYSIS_DIR / "ktypes_sim.csv"
MODIFICATIONS_CSV = ANALYSIS_DIR / "ktypes_modifications.csv"
PLOTS_DIR = ANALYSIS_DIR / "plots"
STRUCTURE_SUBDIR = "plots/ktypes_structures"

TEST_SAMPLE_PATH_KTYPES = Path(__file__).resolve().parent / "sample/ktypes.csv"
TEST_SAMPLE_PATH_MOD = Path(__file__).resolve().parent / "sample/ktypes_modifications.csv"
TEST_SAMPLE_PATH_SIM = Path(__file__).resolve().parent / "sample/ktypes_sim.csv"


def main() -> None:

    ### BUILD TABLES
    tables_api = KTypeTablesAPI(input_xlsx=INPUT_XLSX, output_dir=ANALYSIS_DIR)
    processed_df = tables_api.build_processed_table()
    processed_csv_path = tables_api.export_processed_table(processed_df)

    modifications_api = KTypeModificationsAPI(input_xlsx=INPUT_XLSX, output_dir=ANALYSIS_DIR)
    modifications_df = modifications_api.build_modifications_table()
    modifications_api.export_modifications_table(filename=MODIFICATIONS_CSV.name)

    similarity_api = KTypeSimilarityAPI(input_xlsx=INPUT_XLSX, output_dir=ANALYSIS_DIR)
    similarity_df = similarity_api.build_similarity_table(processed_df)
    similarity_csv_path = similarity_api.export_similarity_table(
        df=similarity_df, filename=SIMILARITY_CSV.name
    )

    #### SAMPLES
    if not processed_df.empty:
        sample_size = min(10, len(processed_df))
        TEST_SAMPLE_PATH_KTYPES.parent.mkdir(parents=True, exist_ok=True)
        processed_df.sample(n=sample_size, random_state=0).to_csv(
            TEST_SAMPLE_PATH_KTYPES, index=False
        )

    if not modifications_df.empty:
        sample_size = min(10, len(modifications_df))
        TEST_SAMPLE_PATH_MOD.parent.mkdir(parents=True, exist_ok=True)
        modifications_df.sample(n=sample_size, random_state=0).to_csv(
            TEST_SAMPLE_PATH_MOD, index=False
        )

 
    if not similarity_df.empty:
        sample_size = min(10, len(similarity_df))
        TEST_SAMPLE_PATH_SIM.parent.mkdir(parents=True, exist_ok=True)
        similarity_df.sample(n=sample_size, random_state=0).to_csv(
            TEST_SAMPLE_PATH_SIM, index=False
        )

    #### DRAW K-TYPE STRUCTURES
    drawer = KTypeStructureDrawer(
        input_xlsx=INPUT_XLSX, output_dir=ANALYSIS_DIR, plots_subdir=STRUCTURE_SUBDIR
    )
    drawer.draw_all(processed_df, test=True)


    #### PLOTS
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    plot_api = KTypePlotAPI(
        ktypes_csv=processed_csv_path,
        similarity_csv=similarity_csv_path,
        modifications_csv=MODIFICATIONS_CSV,
    )

    plot_api.plot_cumulative_structures(
        output_path=PLOTS_DIR / "year.png", show=False
    )

    plot_api.plot_similarity_histogram(
        output_path=PLOTS_DIR / "score_hist.png", column="weighted_total", bins=40, show=False
        )

    plot_api.plot_modification_and_monosaccharide_distribution(
        output_path=PLOTS_DIR / "mono_prevalence.png",
        show=False,
    )
    plot_api.export_similarity_network(
        output_path=PLOTS_DIR / "ktypes_network_THR90.csv", min_weighted_core=0.9
    )


if __name__ == "__main__":
    main()
