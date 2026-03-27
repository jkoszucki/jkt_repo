"""
Minimal standalone PyMOL rendering script.

Requires ONLY PyMOL + stdlib — no yaml, pandas, or matplotlib.

Called as a subprocess from structure_render.py (running in the jkoszucki env):

    conda run -n pymol python structure_render_pymol.py \
        --jobs '[{"cif": "/path/KL1_PC1817.cif", "out": "/path/KL1_PC1817.png"}, ...]' \
        --domain-colors '[["1-350", "#c9a227"]]' \
        --sgnh-color "#c9a227" \
        --dpi 300
"""

from __future__ import annotations

import argparse
import json
import os
import sys

_GRAY_RGB  = [207/255, 207/255, 206/255]   # #cfcfce
_PNG_W, _PNG_H = 3000, 2400


def _hex_to_rgb(hex_color: str) -> list[float]:
    h = hex_color.lstrip("#")
    return [int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4)]


# AlphaFold pLDDT colour scheme (thresholded at 90 / 70 / 50)
_PLDDT_COLORS = [
    (90, [0/255,   83/255,  214/255]),   # very high  > 90  dark blue
    (70, [101/255, 203/255, 243/255]),   # confident  > 70  light blue
    (50, [255/255, 219/255,  19/255]),   # low        > 50  yellow
    ( 0, [255/255, 125/255,  69/255]),   # very low  <= 50  orange
]


def _color_by_plddt(obj: str) -> None:
    """Colour cartoon by b-factor (pLDDT) using the standard AlphaFold scheme."""
    from pymol import cmd
    # Apply from lowest to highest so higher thresholds override
    for i, (threshold, rgb) in enumerate(reversed(_PLDDT_COLORS)):
        name = f"_plddt_{i}"
        cmd.set_color(name, rgb)
        if threshold == 0:
            cmd.color(name, obj)
        else:
            cmd.color(name, f"({obj}) and b > {threshold}")


def _render_one(cif_path: str, png_path: str,
                domain_colors: list[list],
                sgnh_color: str, dpi: int,
                color_by_plddt: bool = False) -> None:
    from pymol import cmd

    os.makedirs(os.path.dirname(png_path), exist_ok=True)

    obj = os.path.splitext(os.path.basename(cif_path))[0]
    cmd.reinitialize()
    cmd.load(cif_path, obj)
    cmd.hide("everything", obj)
    cmd.show("cartoon", obj)

    if color_by_plddt:
        _color_by_plddt(obj)
    else:
        cmd.set_color("_gray_custom", _GRAY_RGB)
        cmd.color("_gray_custom", obj)

        for idx, entry in enumerate(domain_colors):
            residues, color_hex = entry[0], entry[1]
            resolved = color_hex if color_hex is not None else sgnh_color
            color_name = f"_domain_{idx}"
            cmd.set_color(color_name, _hex_to_rgb(resolved))
            sel = f"_sel_{idx}"
            cmd.select(sel, f"({obj}) and resi {residues}")
            cmd.color(color_name, sel)
            cmd.delete(sel)

    cmd.bg_color("white")
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_opaque_background", 1)
    cmd.viewport(_PNG_W, _PNG_H)
    cmd.orient(obj)
    cmd.center(obj)
    cmd.zoom(obj, buffer=50)
    cmd.png(png_path, width=_PNG_W, height=_PNG_H, dpi=dpi, ray=1)
    cmd.delete(obj)


def main() -> None:
    parser = argparse.ArgumentParser(description="Render CIF structures to PNG using PyMOL.")
    parser.add_argument("--jobs", required=True,
                        help='JSON array of {"cif": "...", "out": "..."} objects')
    parser.add_argument("--domain-colors", required=True,
                        help='JSON array of [residue_selection, hex_or_null] pairs')
    parser.add_argument("--sgnh-color", default="#c9a227",
                        help="Fallback colour for None domain entries")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--color-by-plddt", action="store_true",
                        help="Colour by b-factor (pLDDT) using the AlphaFold scheme")
    args = parser.parse_args()

    jobs = json.loads(args.jobs)
    domain_colors = json.loads(args.domain_colors)
    sgnh_color = args.sgnh_color
    dpi = args.dpi
    color_by_plddt = args.color_by_plddt

    ok, failed = 0, []
    for job in jobs:
        cif_path = job["cif"]
        png_path = job["out"]
        try:
            _render_one(cif_path, png_path, domain_colors, sgnh_color, dpi,
                        color_by_plddt=color_by_plddt)
            print(f"  rendered {os.path.basename(cif_path)}", flush=True)
            ok += 1
        except Exception as exc:
            print(f"  FAILED {os.path.basename(cif_path)}: {exc}", flush=True)
            failed.append(cif_path)

    print(f"\nDone: {ok}/{len(jobs)} rendered.")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
