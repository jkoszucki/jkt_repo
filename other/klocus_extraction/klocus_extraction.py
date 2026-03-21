"""
K-locus extraction from annotated GenBank files.

Phase 1: Run Kaptive v3 on FASTA assemblies to detect K-locus coordinates.
Phase 2: Slice those coordinates from GenBank files, preserving CDS annotations.
"""

import json
import subprocess
import sys
from pathlib import Path

_HELPERS_DIR = Path(__file__).resolve().parents[2] / "scripts" / "helpers"

from Bio import SeqIO


class KLocusExtractor:
    def __init__(self, genbank_dir: Path, fasta_dir: Path, kloci_dir: Path, jsonl_path: Path):
        self.genbank_dir = genbank_dir
        self.fasta_dir = fasta_dir
        self.kloci_dir = kloci_dir
        self.jsonl_path = jsonl_path

    def run_kaptive(self):
        fasta_files = sorted(self.fasta_dir.glob("*.fasta"))
        if not fasta_files:
            print(f"ERROR: No FASTA files found in {self.fasta_dir}", file=sys.stderr)
            sys.exit(1)
        self.jsonl_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"Running Kaptive on {len(fasta_files)} assemblies...")
        cmd = [
            "conda", "run", "-n", "kaptive",
            "kaptive", "assembly", "kpsc_k",
            *[str(f) for f in fasta_files],
            "-j", str(self.jsonl_path),
            "-t", "0",
        ]
        result = subprocess.run(cmd, check=True)
        print(f"Kaptive finished (exit {result.returncode}). Results: {self.jsonl_path}")

    def extract(self):
        if not self.jsonl_path.exists():
            print(f"ERROR: {self.jsonl_path} not found. Run run_kaptive() first.", file=sys.stderr)
            sys.exit(1)

        self.kloci_dir.mkdir(parents=True, exist_ok=True)

        with open(self.jsonl_path) as fh:
            lines = [line.strip() for line in fh if line.strip()]

        total = len(lines)
        already_done = sum(1 for l in lines if (self.kloci_dir / f"{json.loads(l)['sample_name']}.gb").exists())
        if already_done:
            print(f"Resuming — {already_done}/{total} already done, skipping.")
        print(f"Extracting K-loci for {total} samples...")
        written = empty = errors = skipped = 0

        for i, line in enumerate(lines, 1):
            data = json.loads(line)
            sample = data["sample_name"]
            best_match = data.get("best_match", "unknown")
            pieces = data.get("pieces", [])
            out_path = self.kloci_dir / f"{sample}.gb"

            if out_path.exists():
                skipped += 1
                print(f"\r  {i}/{total}  {sample:<20} (skip)", end="", flush=True)
                continue

            print(f"\r  {i}/{total}  {sample:<20}      ", end="", flush=True)

            if not pieces:
                out_path.touch()
                empty += 1
                continue

            gb_path = self.genbank_dir / f"{sample}.gb"
            if not gb_path.exists():
                print(f"  WARNING: GenBank not found for {sample}", file=sys.stderr)
                errors += 1
                continue

            try:
                gb_dict = SeqIO.to_dict(SeqIO.parse(gb_path, "genbank"))
                regions = []
                for piece in pieces:
                    contig_id = piece["id"]
                    start = int(piece["start"])
                    end = int(piece["end"])
                    if contig_id not in gb_dict:
                        print(f"  WARNING: contig {contig_id} not found in {sample}.gb", file=sys.stderr)
                        continue
                    region = gb_dict[contig_id][start:end]
                    region.id = f"{sample}_{contig_id}_{start}_{end}"
                    region.name = sample[:16]
                    region.description = f"K-locus region {contig_id}:{start}-{end}"
                    regions.append(region)

                if not regions:
                    out_path.touch()
                    empty += 1
                    continue

                combined = regions[0]
                for r in regions[1:]:
                    combined = combined + r
                combined.id = sample
                combined.name = f"{sample}_{best_match}"
                combined.description = f"{sample} K-locus={best_match}"

                with open(out_path, "w") as out_fh:
                    SeqIO.write(combined, out_fh, "genbank")
                written += 1

            except Exception as exc:
                print(f"  ERROR processing {sample}: {exc}", file=sys.stderr)
                errors += 1

        print()  # newline after progress line
        print(f"Done. Written: {written}, Empty: {empty}, Skipped: {skipped}, Errors: {errors}")


if __name__ == "__main__":
    sys.path.insert(0, str(_HELPERS_DIR))
    from config import Config

    cfg = Config()
    extractor = KLocusExtractor(
        genbank_dir = cfg.input_dir / "gwas/1_BACTERIA/1_GENBANK_GENOMES",
        fasta_dir   = cfg.input_dir / "gwas/1_BACTERIA/2_FASTA_GENOMES_NT",
        kloci_dir   = cfg.input_dir / "gwas/4_K_LOCI",
        jsonl_path  = cfg.output_dir / "preparation" / "kaptive_results.jsonl",
    )
    extractor.run_kaptive()
    extractor.extract()
