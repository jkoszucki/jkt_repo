"""API for accessing K-type modification metadata."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from ktypes_process import KTypeTablesAPI


class KTypeModificationsAPI(KTypeTablesAPI):
    """Expose helper methods focused on the modifications worksheet."""

    MODIFICATIONS_FILENAME = "ktypes_modifications.csv"

    def build_modifications_table(self) -> pd.DataFrame:
        """Return the raw modifications table."""
        return self.get_modifications()

    def export_modifications_table(
        self, df: Optional[pd.DataFrame] = None, filename: Optional[str] = None
    ) -> Path:
        table = df if df is not None else self.build_modifications_table()
        return self.save_dataframe(table, filename or self.MODIFICATIONS_FILENAME)
