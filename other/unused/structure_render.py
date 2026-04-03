"""
Render AF3 CIF structures to PNG using PyMOL.

Reads CIF files from:
    per_locus_dir / {locus} / protein / {locus}_{pc}.cif

Writes PNGs to:
    plots_dir / per-locus / {locus} / {locus}_{pc}.png

Domain colouring is controlled by the domain_colors argument — a list of
(residue_selection, color) tuples passed from figures/main.py:
    residue_selection — PyMOL resi syntax, applied to all chains
    color             — hex string, or None → uses style.sgnh_domain_color
Residues not listed keep the default gray (#cfcfce).

PyMOL is not available in the jkoszucki env, so rendering is delegated to
structure_render_pymol.py via a subprocess call into the pymol conda env.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def render_all_structures(
    per_locus_dir: Path,
    plots_dir: Path,
    domain_colors: list[tuple[str, str | None]],
    style=None,
) -> None:
    """
    Render every CIF in per_locus_dir to a PNG in plots_dir/per-locus/.

    Args:
        per_locus_dir:  cfg.output_dir / "chapter2" / "per-locus"
        plots_dir:      scripts/figures/chapter2/plots
        domain_colors:  [(residue_selection, hex_or_None), ...]
        style:          cfg.style
    """
    sgnh_color = getattr(style, "sgnh_domain_color", "#c9a227")
    dpi        = getattr(style, "dpi", 300)

    cif_files = sorted(per_locus_dir.rglob("*.cif"))
    if not cif_files:
        print("  No CIF files found — skipping.")
        return

    # Build jobs list; skip CIFs whose PNG already exists (checkpoint)
    jobs, skipped = [], []
    for cif in cif_files:
        # Derive output path from the relative path under per_locus_dir,
        # skipping the 'protein' directory level so predictor PNGs land in
        # plots/per-locus/{locus}/ rather than plots/per-locus/protein/.
        # Examples:
        #   {locus}/protein/{stem}.cif   → per-locus/{locus}/{stem}.png
        #   experimental/{sub}/{stem}.cif → per-locus/experimental/{sub}/{stem}.png
        rel_parts = cif.relative_to(per_locus_dir).parts
        dir_parts = [p for p in rel_parts[:-1] if p != "protein"]
        dest_dir = plots_dir / "per-locus" / Path(*dir_parts)
        png_path = dest_dir / f"{cif.stem}.png"
        if png_path.exists():
            skipped.append(cif.stem)
        else:
            jobs.append({"cif": str(cif), "out": str(png_path)})

    print(f"  Checkpoint: {len(skipped)} already rendered, {len(jobs)} missing.")
    if skipped:
        for s in skipped:
            print(f"    [skip] {s}")

    if not jobs:
        print("  All structures already rendered — nothing to do.")
        return

    # Resolve None → sgnh_color before serialising
    resolved_domain_colors = [
        [resi, (color if color is not None else sgnh_color)]
        for resi, color in domain_colors
    ]

    script_path = Path(__file__).resolve().parent / "structure_render_pymol.py"

    cmd = [
        "conda", "run", "-n", "pymol",
        sys.executable if False else "python",   # use env's python
        str(script_path),
        "--jobs", json.dumps(jobs),
        "--domain-colors", json.dumps(resolved_domain_colors),
        "--sgnh-color", sgnh_color,
        "--dpi", str(dpi),
    ]

    print(f"  Rendering {len(jobs)} structure(s) via pymol env …")
    result = subprocess.run(cmd, text=True, capture_output=False)
    if result.returncode != 0:
        print(f"  [warn] structure_render_pymol.py exited with code {result.returncode}")
    else:
        print(f"  Structure PNGs → {plots_dir / 'per-locus'}")


def render_af3_no_ecod(
    af3_folds_dir: Path,
    pcs: list[str],
    plots_dir: Path,
    style=None,
) -> None:
    """
    Render model_0 CIF from AF3 raw predictions for each PC in *pcs*.

    Args:
        af3_folds_dir:  cfg.output_dir / "chapter2" / "raw_af3" / "folds_YYYY_MM_DD_HH_MM"
        pcs:            list of PC identifiers, e.g. ["PC0232", "PC0294", ...]
        plots_dir:      scripts/figures/chapter2/plots
        style:          cfg.style
    """
    dpi = getattr(style, "dpi", 300)

    out_dir = plots_dir / "af3_no_ecod"
    out_dir.mkdir(parents=True, exist_ok=True)

    jobs, skipped, missing = [], [], []
    for pc in pcs:
        pc_lower = pc.lower()
        cif = af3_folds_dir / pc_lower / f"fold_{pc_lower}_model_0.cif"
        if not cif.exists():
            missing.append(pc)
            continue
        png = out_dir / f"{pc}.png"
        if png.exists():
            skipped.append(pc)
        else:
            jobs.append({"cif": str(cif), "out": str(png)})

    if missing:
        for pc in missing:
            print(f"  [warn] model_0 CIF not found for {pc} — skipping")
    if skipped:
        print(f"  Checkpoint: {len(skipped)} already rendered — {', '.join(skipped)}")
    if not jobs:
        print("  All AF3 no-ECOD structures already rendered — nothing to do.")
        return

    script_path = Path(__file__).resolve().parent / "structure_render_pymol.py"
    cmd = [
        "conda", "run", "-n", "pymol", "python",
        str(script_path),
        "--jobs", json.dumps(jobs),
        "--domain-colors", "[]",
        "--dpi", str(dpi),
    ]

    print(f"  Rendering {len(jobs)} AF3 no-ECOD structure(s) via pymol env …")
    result = subprocess.run(cmd, text=True, capture_output=False)
    if result.returncode != 0:
        print(f"  [warn] structure_render_pymol.py exited with code {result.returncode}")
    else:
        print(f"  AF3 no-ECOD PNGs → {out_dir}")


def render_af3_no_ecod_plddt(
    af3_folds_dir: Path,
    pcs: list[str],
    plots_dir: Path,
    style=None,
) -> None:
    """
    Render model_0 CIF coloured by pLDDT (b-factor) for each PC in *pcs*.

    Uses the standard AlphaFold colour scheme:
      >90  dark blue  (very high confidence)
      >70  light blue (confident)
      >50  yellow     (low confidence)
      ≤50  orange     (very low confidence)

    Args:
        af3_folds_dir:  cfg.output_dir / "chapter2" / "raw_af3" / "folds_YYYY_MM_DD_HH_MM"
        pcs:            list of PC identifiers, e.g. ["PC0232", "PC0294", ...]
        plots_dir:      scripts/figures/chapter2/plots
        style:          cfg.style
    """
    dpi = getattr(style, "dpi", 300)

    out_dir = plots_dir / "af3_no_ecod_plddt"
    out_dir.mkdir(parents=True, exist_ok=True)

    jobs, skipped, missing = [], [], []
    for pc in pcs:
        pc_lower = pc.lower()
        cif = af3_folds_dir / pc_lower / f"fold_{pc_lower}_model_0.cif"
        if not cif.exists():
            missing.append(pc)
            continue
        png = out_dir / f"{pc}.png"
        if png.exists():
            skipped.append(pc)
        else:
            jobs.append({"cif": str(cif), "out": str(png)})

    if missing:
        for pc in missing:
            print(f"  [warn] model_0 CIF not found for {pc} — skipping")
    if skipped:
        print(f"  Checkpoint: {len(skipped)} already rendered — {', '.join(skipped)}")
    if not jobs:
        print("  All AF3 pLDDT structures already rendered — nothing to do.")
        return

    script_path = Path(__file__).resolve().parent / "structure_render_pymol.py"
    cmd = [
        "conda", "run", "-n", "pymol", "python",
        str(script_path),
        "--jobs", json.dumps(jobs),
        "--domain-colors", "[]",
        "--dpi", str(dpi),
        "--color-by-plddt",
    ]

    print(f"  Rendering {len(jobs)} AF3 pLDDT structure(s) via pymol env …")
    result = subprocess.run(cmd, text=True, capture_output=False)
    if result.returncode != 0:
        print(f"  [warn] structure_render_pymol.py exited with code {result.returncode}")
    else:
        print(f"  AF3 pLDDT PNGs → {out_dir}")
