"""
K-locus pairwise wGRR similarity — all isolates vs all.

Strategy: single all-vs-all BLASTP (one database + one query), then
compute wGRR from bidirectional best hits + powerneedle identities.

Steps:
1. Merge all per-isolate .faa files (extracted from 4_K_LOCI/*.gb)
2. makeblastdb once
3. blastp all-vs-all (one run)
4. Find BBH (bidirectional best hits, e-value < 0.1)
5. powerneedle on BBH pairs → percent identity
6. wGRR = sum(identity/100) / min(genes1, genes2)
7. Write grr_results.csv
"""

import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

_HELPERS_DIR = Path(__file__).resolve().parents[2] / "scripts" / "helpers"

from Bio import SeqIO

EVALUE  = 0.1
WORKERS = 8


# ---------------------------------------------------------------------------
# .gb → .faa extraction
# ---------------------------------------------------------------------------

def extract_faa(gb_path, faa_path):
    """Extract CDS translations from a GenBank file → multifasta .faa."""
    records = []
    for rec in SeqIO.parse(gb_path, "genbank"):
        for feat in rec.features:
            if feat.type != "CDS" or "translation" not in feat.qualifiers:
                continue
            product = feat.qualifiers.get("product", ["unknown"])[0]
            seq_id  = "{}_{}_{}".format(
                rec.id,
                int(feat.location.start),
                int(feat.location.end),
            )
            aa = feat.qualifiers["translation"][0]
            lines = [f">{seq_id} {product}"] + [aa[i:i+60] for i in range(0, len(aa), 60)]
            records.append("\n".join(lines))
    with open(faa_path, "w") as f:
        f.write("\n".join(records) + "\n")
    return len(records)


def prepare_faa_files(k_loci_dir, faa_dir):
    """Extract .faa for every non-empty .gb; return {sample: faa_path}."""
    faa_dir.mkdir(parents=True, exist_ok=True)
    faa_paths = {}
    gb_files = sorted(p for p in k_loci_dir.glob("*.gb") if p.stat().st_size > 0)
    print(f"  {len(gb_files)} non-empty .gb files")
    for gb in gb_files:
        faa = faa_dir / (gb.stem + ".faa")
        if not faa.exists():
            extract_faa(gb, faa)
        if faa.stat().st_size > 0:
            faa_paths[gb.stem] = faa
    print(f"  {len(faa_paths)} samples with protein sequences")
    return faa_paths


# ---------------------------------------------------------------------------
# Step 1 – merge all .faa, tag protein IDs as SAMPLE__PROTID
# ---------------------------------------------------------------------------

def merge_faa(faa_dir, merged_path):
    n_samples = n_proteins = 0
    with open(merged_path, "w") as out:
        for faa in sorted(faa_dir.glob("*.faa")):
            sample = faa.stem
            for line in faa.read_text().splitlines():
                if line.startswith(">"):
                    prot_id = line[1:].split()[0]
                    out.write(f">{sample}__{prot_id}\n")
                    n_proteins += 1
                else:
                    out.write(line + "\n")
            n_samples += 1
    print(f"  {n_samples} samples, {n_proteins:,} proteins → {merged_path.name}")


# ---------------------------------------------------------------------------
# Step 2 – makeblastdb
# ---------------------------------------------------------------------------

def make_blast_db(merged_faa, db_prefix):
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        f"makeblastdb -in {merged_faa} -dbtype prot -out {db_prefix}",
        shell=True, check=True,
    )


# ---------------------------------------------------------------------------
# Step 3 – all-vs-all blastp
# ---------------------------------------------------------------------------

def run_blast(merged_faa, db_prefix, out_tsv, threads=WORKERS):
    subprocess.run(
        f"blastp -query {merged_faa} -db {db_prefix} "
        f"-outfmt '7 qseqid sseqid evalue pident' "
        f"-evalue {EVALUE} -num_threads {threads} "
        f"-out {out_tsv}",
        shell=True, check=True,
    )


# ---------------------------------------------------------------------------
# Step 4 – parse blast → BBH
# ---------------------------------------------------------------------------

def parse_blast_bbh(blast_tsv):
    """Return {(sampleA, sampleB): [(protA, protB)]} with A < B."""
    best = {}   # (q_sample, q_prot, s_sample) -> (s_prot, evalue)
    print("  Parsing BLAST output…")
    with open(blast_tsv) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            qid, sid = parts[0], parts[1]
            evalue   = float(parts[2])
            q_sample, q_prot = qid.split("__", 1)
            s_sample, s_prot = sid.split("__", 1)
            if q_sample == s_sample:
                continue
            key = (q_sample, q_prot, s_sample)
            if key not in best or evalue < best[key][1]:
                best[key] = (s_prot, evalue)
    print(f"  {len(best):,} best inter-sample hits")

    bbh = defaultdict(list)
    for (q_sample, q_prot, s_sample), (s_prot, _) in best.items():
        if best.get((s_sample, s_prot, q_sample), (None,))[0] == q_prot:
            a = (q_sample, q_prot) if q_sample <= s_sample else (s_sample, s_prot)
            b = (s_sample, s_prot) if q_sample <= s_sample else (q_sample, q_prot)
            bbh[(a[0], b[0])].append((a[1], b[1]))
    print(f"  {len(bbh):,} pairs have ≥1 BBH")
    return bbh


# ---------------------------------------------------------------------------
# Step 5 – powerneedle per pair
# ---------------------------------------------------------------------------

def _write_pair_faa(sample_a, sample_b, prots_a, prots_b, faa_dir, out_path):
    needed = {sample_a: set(prots_a), sample_b: set(prots_b)}
    with open(out_path, "w") as out:
        for sample in (sample_a, sample_b):
            emit = False
            current_needed = needed[sample]
            for line in (faa_dir / f"{sample}.faa").read_text().splitlines():
                if line.startswith(">"):
                    emit = line[1:].split()[0] in current_needed
                if emit:
                    out.write(line + "\n")


def _run_needle_pair(args):
    sample_a, sample_b, bbh_pairs, faa_dir, tmp_dir = args
    tmp_dir.mkdir(exist_ok=True)
    base   = tmp_dir / f"{sample_a}-{sample_b}"
    faa_p  = Path(str(base) + ".faa")
    bbh_p  = Path(str(base) + ".bbh")
    ndl_p  = Path(str(base) + ".needle")
    alg_p  = Path(str(base) + ".needlealg")

    _write_pair_faa(sample_a, sample_b,
                    [p[0] for p in bbh_pairs],
                    [p[1] for p in bbh_pairs],
                    faa_dir, faa_p)
    bbh_p.write_text("".join(f"{a}\t{b}\n" for a, b in bbh_pairs))

    subprocess.run(
        f"powerneedle -gapopen 10 -gapextend 0.5 -brief "
        f"-pairs {bbh_p} -identities {ndl_p} -alignment {alg_p} {faa_p}",
        shell=True, capture_output=True,
    )

    identities = []
    if ndl_p.exists():
        for line in ndl_p.read_text().splitlines():
            if not line.startswith("Sequence1"):
                try:
                    identities.append(float(line.split()[2]))
                except (IndexError, ValueError):
                    pass

    for p in (faa_p, bbh_p, ndl_p, alg_p):
        p.unlink(missing_ok=True)

    return sample_a, sample_b, identities


# ---------------------------------------------------------------------------
# Step 6 – assemble output CSV
# ---------------------------------------------------------------------------

def _gene_counts(faa_dir):
    return {
        faa.stem: sum(1 for l in faa.read_text().splitlines() if l.startswith(">"))
        for faa in faa_dir.glob("*.faa")
    }


def _compute_write(bbh, faa_dir, tmp_dir, out_csv):
    gene_counts = _gene_counts(faa_dir)
    total = len(bbh)
    print(f"  Running powerneedle for {total:,} pairs ({WORKERS} workers)…")

    rows = []
    done = 0
    args_list = [(a, b, pairs, faa_dir, tmp_dir) for (a, b), pairs in sorted(bbh.items())]

    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures = {pool.submit(_run_needle_pair, arg): arg[:2] for arg in args_list}
        for fut in as_completed(futures):
            sample_a, sample_b, identities = fut.result()
            g1 = gene_counts.get(sample_a, 0)
            g2 = gene_counts.get(sample_b, 0)
            grr_sum = sum(i / 100 for i in identities)
            wgrr    = grr_sum / min(g1, g2) if min(g1, g2) > 0 else 0
            min35   = sum(1 for i in identities if i >= 35)
            min80   = sum(1 for i in identities if i >= 80)
            rows.append(
                f"{sample_a}\t{sample_b}\t{wgrr:.6f}\t{grr_sum:.4f}"
                f"\t{g1}\t{g2}\t{len(bbh[(sample_a,sample_b)])}\t{min35}\t{min80}"
            )
            done += 1
            if done % 50_000 == 0 or done == total:
                print(f"    {done:,}/{total:,} done", flush=True)

    header = "sample1\tsample2\twgrr\tgrr_sum\tgenes1\tgenes2\tbbh\tmin35\tmin80"
    out_csv.write_text(header + "\n" + "\n".join(rows) + "\n")
    print(f"  {len(rows):,} rows → {out_csv.name}")


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def compute_klocus_grr(k_loci_dir, output_dir):
    """
    Compute all-vs-all wGRR for K-loci.

    Args:
        k_loci_dir:  Path to input_dir/gwas/4_K_LOCI/
        output_dir:  Path to cfg.output_dir / "chapter2"

    Writes:
        output_dir / grr_results.csv
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir    = output_dir / "grr_workdir"
    faa_dir     = work_dir / "faa"
    merged_faa  = work_dir / "all_proteins.faa"
    blast_db    = work_dir / "blastdb" / "kloci"
    blast_out   = work_dir / "blast_all_vs_all.tsv"
    tmp_dir     = work_dir / "tmp_needle"
    out_csv     = output_dir / "grr_results.csv"

    print("Step 1: Extracting protein sequences…")
    prepare_faa_files(k_loci_dir, faa_dir)

    print("Step 2: Merging .faa files…")
    if not merged_faa.exists():
        merge_faa(faa_dir, merged_faa)
    else:
        print(f"  {merged_faa.name} already exists, skipping")

    print("Step 3: Building BLAST database…")
    db_check = blast_db.parent / (blast_db.name + ".pin")
    if not db_check.exists():
        make_blast_db(merged_faa, blast_db)
    else:
        print("  DB already exists, skipping")

    print("Step 4: Running all-vs-all BLASTP…")
    if not blast_out.exists():
        run_blast(merged_faa, blast_db, blast_out)
    else:
        print(f"  {blast_out.name} already exists, skipping")

    print("Step 5: Finding BBH and computing wGRR…")
    bbh = parse_blast_bbh(blast_out)
    _compute_write(bbh, faa_dir, tmp_dir, out_csv)

    print("Done.")


if __name__ == "__main__":
    sys.path.insert(0, str(_HELPERS_DIR))
    from config import Config

    cfg = Config()
    compute_klocus_grr(
        k_loci_dir = cfg.input_dir / "gwas/4_K_LOCI",
        output_dir = cfg.output_dir / "chapter2",
    )
