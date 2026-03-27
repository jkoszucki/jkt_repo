"""
K-locus pairwise wGRR similarity — all isolates vs all.

Follows GRRpair methodology exactly (Kupczok et al. 2022 / de Sousa et al.
2021): same BLAST format, same BBH extraction, same powerneedle call, same
wGRR formula.

GRRpair repo: /Users/januszkoszucki/Projects/code/GRRpair

Scalability note: GRRpair runs BLAST per pair, which is infeasible for ~3,911
isolates (~7.6 M pairs). Here we run a single all-vs-all BLASTP, then extract
per-pair BBH from it — the BBH logic is identical to GRRpair's extract_best.

Steps:
1. Extract .faa per non-empty K-locus .gb file
2. Merge all .faa, tagging protein IDs as SAMPLE__PROTID
3. makeblastdb once
4. blastp all-vs-all (outfmt 7, same as GRRpair)
5. Extract BBH per pair (GRRpair's extract_best logic)
6. powerneedle on BBH pairs (same flags as GRRpair)
7. wGRR = sum(identity/100) / min(genes1, genes2)  (GRRpair formula)
8. Write grr_results.csv
"""

import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO

EVALUE  = 0.1
WORKERS = 8


# ---------------------------------------------------------------------------
# Step 1 – .gb → .faa extraction
# ---------------------------------------------------------------------------

def extract_faa(gb_path, faa_path):
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
            aa    = feat.qualifiers["translation"][0]
            lines = [f">{seq_id} {product}"] + [aa[i:i+60] for i in range(0, len(aa), 60)]
            records.append("\n".join(lines))
    with open(faa_path, "w") as f:
        f.write("\n".join(records) + "\n")
    return len(records)


def prepare_faa_files(k_loci_dir, faa_dir):
    faa_dir.mkdir(parents=True, exist_ok=True)
    faa_paths = {}
    gb_files  = sorted(p for p in k_loci_dir.glob("*.gb") if p.stat().st_size > 0)
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
# Step 2 – merge all .faa, tagging IDs as SAMPLE__PROTID
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
# Step 3 – makeblastdb
# ---------------------------------------------------------------------------

def make_blast_db(merged_faa, db_prefix):
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        f"makeblastdb -in {merged_faa} -dbtype prot -out {db_prefix}",
        shell=True, check=True,
    )


# ---------------------------------------------------------------------------
# Step 4 – all-vs-all blastp (outfmt 7 — same as GRRpair)
# ---------------------------------------------------------------------------

def run_blast(merged_faa, db_prefix, out_tsv, threads=WORKERS):
    subprocess.run(
        f"blastp -query {merged_faa} -db {db_prefix} "
        f"-outfmt 7 -evalue {EVALUE} -num_threads {threads} "
        f"-out {out_tsv}",
        shell=True, check=True,
    )


# ---------------------------------------------------------------------------
# Step 5 – BBH extraction (GRRpair's extract_best logic, adapted for
#           all-vs-all: track best hit per (q_sample, q_prot, s_sample))
# ---------------------------------------------------------------------------

def parse_blast_bbh(blast_tsv):
    """
    Mirrors GRRpair's extract_best: take the first (best e-value) hit for
    each query protein within each target sample, filter by evalue < 0.1.
    Returns {(sampleA, sampleB): [(protA, protB)]} with A <= B.
    """
    # best[(q_sample, q_prot, s_sample)] = s_prot  (first = best e-value hit)
    best = {}
    print("  Parsing BLAST output…")
    with open(blast_tsv) as f:
        for line in f:
            if line.startswith("#"):
                continue
            spl = line.split()
            if len(spl) < 11:
                continue
            qid    = spl[0]
            sid    = spl[1]
            evalue = float(spl[10])          # column 10 in outfmt 7
            if evalue >= EVALUE:
                continue
            q_sample, q_prot = qid.split("__", 1)
            s_sample, s_prot = sid.split("__", 1)
            if q_sample == s_sample:
                continue
            key = (q_sample, q_prot, s_sample)
            if key not in best:              # first occurrence = best e-value
                best[key] = s_prot

    print(f"  {len(best):,} best inter-sample hits")

    bbh = defaultdict(list)
    seen = set()
    for (q_sample, q_prot, s_sample), s_prot in best.items():
        if best.get((s_sample, s_prot, q_sample)) == q_prot:
            a_sample = min(q_sample, s_sample)
            b_sample = max(q_sample, s_sample)
            a_prot   = q_prot if q_sample == a_sample else s_prot
            b_prot   = s_prot if q_sample == a_sample else q_prot
            pair_key = (a_sample, a_prot, b_sample, b_prot)
            if pair_key not in seen:
                seen.add(pair_key)
                bbh[(a_sample, b_sample)].append((a_prot, b_prot))

    print(f"  {len(bbh):,} pairs have ≥1 BBH")
    return bbh


# ---------------------------------------------------------------------------
# Step 6+7 – powerneedle + wGRR per pair (GRRpair's exact procedure)
# ---------------------------------------------------------------------------

def _run_pair(args):
    """
    For one (sampleA, sampleB) pair:
      - cat faa_a faa_b  →  pair.faa          (GRRpair: cat genome1.faa genome2.faa)
      - write .bbh file                        (GRRpair: outf.write(...))
      - powerneedle -gapopen 10 -gapextend 0.5 -brief ...  (identical flags)
      - parse .needle identities               (GRRpair: line.split()[2])
      - wGRR = sum(ident/100) / min(g1, g2)   (GRRpair formula)
    """
    sample_a, sample_b, bbh_pairs, faa_dir, tmp_dir = args
    tmp_dir.mkdir(exist_ok=True)

    faa_a = faa_dir / f"{sample_a}.faa"
    faa_b = faa_dir / f"{sample_b}.faa"
    base  = tmp_dir / f"{sample_a}-{sample_b}"
    faa_p = Path(str(base) + ".faa")
    bbh_p = Path(str(base) + ".bbh")
    ndl_p = Path(str(base) + ".needle")
    alg_p = Path(str(base) + ".needlealg")

    # cat genome1.faa genome2.faa > pair.faa  (GRRpair line 58)
    subprocess.run(f"cat {faa_a} {faa_b} > {faa_p}", shell=True, check=True)

    # write .bbh  (GRRpair lines 53-56)
    bbh_p.write_text("".join(f"{a}\t{b}\n" for a, b in bbh_pairs))

    # powerneedle  (GRRpair line 59)
    subprocess.run(
        f"powerneedle -gapopen 10 -gapextend 0.5 -brief "
        f"-pairs {bbh_p} -identities {ndl_p} -alignment {alg_p} {faa_p}",
        shell=True, capture_output=True,
    )

    # gene counts  (GRRpair lines 64-65)
    g1 = sum(1 for l in faa_a.read_text().splitlines() if l.startswith(">"))
    g2 = sum(1 for l in faa_b.read_text().splitlines() if l.startswith(">"))

    # parse identities  (GRRpair lines 70-75)
    grr = min35 = min80 = 0
    if ndl_p.exists():
        for line in ndl_p.read_text().splitlines():
            if not line.startswith("Sequence1"):
                try:
                    ident = float(line.split()[2])
                    grr  += ident / 100
                    if ident >= 35: min35 += 1
                    if ident >= 80: min80 += 1
                except (IndexError, ValueError):
                    pass

    wgrr = grr / min(g1, g2) if min(g1, g2) > 0 else 0

    for p in (faa_p, bbh_p, ndl_p, alg_p):
        p.unlink(missing_ok=True)

    return sample_a, sample_b, wgrr, grr, g1, g2, len(bbh_pairs), min35, min80


def compute_write_grr(bbh, faa_dir, tmp_dir, out_csv):
    total = len(bbh)
    print(f"  Running powerneedle for {total:,} pairs ({WORKERS} workers)…")

    rows = []
    done = 0
    args_list = [
        (a, b, pairs, faa_dir, tmp_dir)
        for (a, b), pairs in sorted(bbh.items())
    ]

    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures = {pool.submit(_run_pair, arg): arg[:2] for arg in args_list}
        for fut in as_completed(futures):
            sample_a, sample_b, wgrr, grr, g1, g2, n_bbh, min35, min80 = fut.result()
            rows.append(
                f"{sample_a}\t{sample_b}\t{wgrr:.6f}\t{grr:.4f}"
                f"\t{g1}\t{g2}\t{n_bbh}\t{min35}\t{min80}"
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

def compute_klocus_grr(k_loci_dir, work_dir, grr_output_path):
    """
    Compute all-vs-all wGRR for K-loci.

    Args:
        k_loci_dir:       Path to input_dir/gwas/4_K_LOCI/
        work_dir:         Scratch directory for BLAST DB, temp needle files, etc.
                          (large; keep outside the repository)
        grr_output_path:  Destination for grr_results.csv
                          (written to input_dir/gwas/ so downstream scripts can
                          read it as raw input — same convention as 4_K_LOCI/)

    Writes:
        grr_output_path   (e.g. input_dir/gwas/grr_results.csv)
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    faa_dir    = work_dir / "faa"
    merged_faa = work_dir / "all_proteins.faa"
    blast_db   = work_dir / "blastdb" / "kloci"
    blast_out  = work_dir / "blast_all_vs_all.tsv"
    tmp_dir    = work_dir / "tmp_needle"
    out_csv    = grr_output_path

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

    print("Step 5: Extracting BBH (GRRpair logic)…")
    bbh = parse_blast_bbh(blast_out)

    print("Step 6: powerneedle + wGRR (GRRpair formula)…")
    compute_write_grr(bbh, faa_dir, tmp_dir, out_csv)

    print("Done.")


if __name__ == "__main__":
    _HELPERS_DIR = Path(__file__).resolve().parents[2] / "scripts" / "helpers"
    sys.path.insert(0, str(_HELPERS_DIR))
    from config import Config

    cfg = Config()
    compute_klocus_grr(
        k_loci_dir      = cfg.input_dir / "gwas/4_K_LOCI",
        work_dir        = cfg.output_dir / "grr_workdir",
        grr_output_path = cfg.input_dir / "gwas/grr_results.csv",
    )
