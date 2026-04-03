"""
Step 5: Mirror GWAS data for SSLBH-positive GWAS PCs.

Creates symlinks in acetyl_dir/acetyl-gwas/{locus}/{clustering_level}/{PC}/
pointing to the corresponding gwas-data/ hierarchy directory.

Symlink target: gwas_root/{ecod_folder}/{locus}/{clustering_level}/{PC}/
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def mirror_gwas_data(
    acetyl_gwas_tsv: Path,
    gwas_root: Path,
    acetyl_dir: Path,
) -> None:
    """
    Create symlinks for each SSLBH-positive GWAS predictor.

    acetyl_dir/acetyl-gwas/{locus}/{clustering_level}/{PC}/
        → gwas_root/{ecod_folder}/{locus}/{clustering_level}/{PC}/

    Skips if symlink already exists and is valid.
    """
    df = pd.read_csv(acetyl_gwas_tsv, sep="\t")
    mirror_root = acetyl_dir / "acetyl-gwas"
    mirror_root.mkdir(parents=True, exist_ok=True)

    ok = skipped = broken = 0

    for _, row in df.iterrows():
        locus = row["locus"]
        cl    = row["clustering_level"]
        pc    = row["PC"]
        ecod_folder = row["ecod_folder"]

        target = gwas_root / ecod_folder / locus / cl / pc
        link   = mirror_root / locus / cl / pc

        if not target.exists():
            print(f"  [warn] target not found, skipping: {target}")
            broken += 1
            continue

        link.parent.mkdir(parents=True, exist_ok=True)

        if link.is_symlink():
            if link.resolve() == target.resolve():
                skipped += 1
                continue
            link.unlink()  # stale symlink — replace

        link.symlink_to(target)
        ok += 1

    total = len(df)
    print(
        f"  Symlinks: {ok} created, {skipped} already valid, "
        f"{broken} skipped (target missing) — {total} rows total"
    )
