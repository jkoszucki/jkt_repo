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
        """Return the modifications table with structure_id joined from the ktypes sheet."""
        mods = self.get_modifications()
        raw = self.get_sheet(self.KTYPE_SHEET)
        if "structure_id" in raw.columns and "structure_id" not in mods.columns:
            k_col = next(
                (c for c in raw.columns if str(c).strip().lower().replace("-", "_") in {"k_type", "ktype"}),
                None,
            )
            mod_k_col = next(
                (c for c in mods.columns if str(c).strip().lower().replace("-", "_") in {"k_type", "ktype"}),
                None,
            )
            if k_col and mod_k_col:
                id_map = dict(zip(
                    raw[k_col].astype(str).str.strip(),
                    raw["structure_id"].astype(str).str.strip(),
                ))
                mods.insert(0, "structure_id", mods[mod_k_col].astype(str).str.strip().map(id_map))
        return mods

    def export_modifications_table(
        self, df: Optional[pd.DataFrame] = None, filename: Optional[str] = None
    ) -> Path:
        table = df if df is not None else self.build_modifications_table()
        return self.save_dataframe(table, filename or self.MODIFICATIONS_FILENAME)
