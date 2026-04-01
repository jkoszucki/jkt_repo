from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from ktypes_modifications import KTypeModificationsAPI
from ktypes_process import KTypeTablesAPI
from ktypes_similarity import KTypeSimilarityAPI

SAMPLE_DIR = Path(__file__).resolve().parent / "sample"


def main() -> None:
    cfg = Config()
    input_xlsx = cfg.input_dir / "cps.xlsx"
    chapter3_dir = cfg.output_dir / "chapter3"

    ### BUILD TABLES
    tables_api = KTypeTablesAPI(input_xlsx=input_xlsx, output_dir=chapter3_dir)
    processed_df = tables_api.build_processed_table()
    tables_api.export_processed_table(processed_df)

    modifications_api = KTypeModificationsAPI(input_xlsx=input_xlsx, output_dir=chapter3_dir)
    modifications_df = modifications_api.build_modifications_table()
    modifications_api.export_modifications_table()

    similarity_api = KTypeSimilarityAPI(input_xlsx=input_xlsx, output_dir=chapter3_dir)
    similarity_df = similarity_api.build_similarity_table(processed_df)
    similarity_api.export_similarity_table(df=similarity_df)

    ### SAMPLES
    if not processed_df.empty:
        sample_size = min(10, len(processed_df))
        SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
        processed_df.sample(n=sample_size, random_state=0).to_csv(
            SAMPLE_DIR / "ktypes.csv", index=False
        )

    if not modifications_df.empty:
        sample_size = min(10, len(modifications_df))
        SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
        modifications_df.sample(n=sample_size, random_state=0).to_csv(
            SAMPLE_DIR / "ktypes_modifications.csv", index=False
        )

    if not similarity_df.empty:
        sample_size = min(10, len(similarity_df))
        SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
        similarity_df.sample(n=sample_size, random_state=0).to_csv(
            SAMPLE_DIR / "ktypes_sim.csv", index=False
        )


if __name__ == "__main__":
    main()
