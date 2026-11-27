#!/usr/bin/env python3

from pathlib import Path
import subprocess
import re
import csv
import os
from itertools import combinations
from statistics import mean

SCRIPT_DIR = Path(__file__).resolve().parent
BASE = SCRIPT_DIR.parent

os.chdir(BASE)

INPUT_FASTA = Path("input/halobacteria_13_score_positive.fasta")
MEMBERSHIP_TSV = Path("tables/cluster_membership_n13_cdhit_0.90.tsv")

PAIRWISE_DIR = Path("pairwise")
PAIR_DIR = PAIRWISE_DIR / "individual_pairs"
TABLES_DIR = Path("tables")
LOGS_DIR = Path("logs")

PAIRWISE_DIR.mkdir(exist_ok=True)
PAIR_DIR.mkdir(exist_ok=True)
TABLES_DIR.mkdir(exist_ok=True)
LOGS_DIR.mkdir(exist_ok=True)

PAIRWISE_TSV = TABLES_DIR / "pairwise_global_identity_n13_cdhit_0.90.tsv"
SUMMARY_TSV = TABLES_DIR / "pairwise_global_identity_summary_n13_cdhit_0.90.tsv"


def read_fasta(path):
    seqs = {}
    current_id = None
    parts = []

    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(parts)
                current_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())

    if current_id is not None:
        seqs[current_id] = "".join(parts)

    return seqs


def read_membership(path):
    clusters = {}

    with path.open() as f:
        next(f)
        for line in f:
            if not line.strip():
                continue
            cluster_id, member_count, members = line.rstrip("\n").split("\t")
            clusters[cluster_id] = members.split(",")

    return clusters


def write_single_fasta(path, acc, seq):
    with path.open("w") as f:
        f.write(f">{acc}\n{seq}\n")


def parse_pair_output(path):
    seq1 = None
    seq2 = None

    re_seq1 = re.compile(r"^#\s*1:\s*(\S+)")
    re_seq2 = re.compile(r"^#\s*2:\s*(\S+)")
    re_identity = re.compile(r"^#\s*Identity:\s*(\d+)/(\d+)\s*\(([\d.]+)%\)")

    with path.open(errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")

            m1 = re_seq1.match(line)
            if m1:
                seq1 = m1.group(1)
                continue

            m2 = re_seq2.match(line)
            if m2:
                seq2 = m2.group(1)
                continue

            mi = re_identity.match(line)
            if mi:
                if seq1 is None or seq2 is None:
                    raise RuntimeError(f"Identity found before sequence IDs in {path}")

                return {
                    "seq1": seq1,
                    "seq2": seq2,
                    "matches": int(mi.group(1)),
                    "aligned_length": int(mi.group(2)),
                    "identity_percent": float(mi.group(3)),
                }

    raise RuntimeError(f"No '# Identity:' field parsed from {path}")


def run_pair(cluster_id, acc1, acc2, seqs):
    safe_name = f"{cluster_id}__{acc1}__{acc2}".replace(".", "_")
    fa1 = PAIR_DIR / f"{safe_name}.A.fasta"
    fa2 = PAIR_DIR / f"{safe_name}.B.fasta"
    out = PAIR_DIR / f"{safe_name}.needle.txt"
    log = LOGS_DIR / f"needleall_{safe_name}.log"

    write_single_fasta(fa1, acc1, seqs[acc1])
    write_single_fasta(fa2, acc2, seqs[acc2])

    cmd = [
        "needleall",
        "-asequence", str(fa1),
        "-bsequence", str(fa2),
        "-gapopen", "10.0",
        "-gapextend", "0.5",
        "-aformat3", "pair",
        "-outfile", str(out),
        "-auto",
    ]

    with log.open("w") as lf:
        subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, check=True)

    record = parse_pair_output(out)

    # Preserve the cluster-defined order, regardless of how EMBOSS reports IDs.
    record["cluster_id"] = cluster_id
    record["seq1"] = acc1
    record["seq2"] = acc2

    return record


def main():
    seqs = read_fasta(INPUT_FASTA)
    clusters = read_membership(MEMBERSHIP_TSV)

    all_rows = []
    summary_rows = []

    for cluster_id, members in clusters.items():
        n = len(members)

        if n == 1:
            summary_rows.append({
                "cluster_id": cluster_id,
                "member_count": n,
                "pair_count": 0,
                "min_identity_percent": "NA",
                "mean_identity_percent": "NA",
                "max_identity_percent": "NA",
                "members": ",".join(members),
            })
            continue

        records = []
        for acc1, acc2 in combinations(members, 2):
            if acc1 not in seqs:
                raise KeyError(f"{acc1} not found in FASTA")
            if acc2 not in seqs:
                raise KeyError(f"{acc2} not found in FASTA")

            records.append(run_pair(cluster_id, acc1, acc2, seqs))

        expected_pairs = n * (n - 1) // 2
        if len(records) != expected_pairs:
            raise RuntimeError(f"{cluster_id}: expected {expected_pairs} pairs, got {len(records)}")

        identities = [r["identity_percent"] for r in records]

        for r in records:
            all_rows.append({
                "cluster_id": r["cluster_id"],
                "seq1": r["seq1"],
                "seq2": r["seq2"],
                "matches": r["matches"],
                "aligned_length": r["aligned_length"],
                "identity_percent": f"{r['identity_percent']:.3f}",
            })

        summary_rows.append({
            "cluster_id": cluster_id,
            "member_count": n,
            "pair_count": len(records),
            "min_identity_percent": f"{min(identities):.3f}",
            "mean_identity_percent": f"{mean(identities):.3f}",
            "max_identity_percent": f"{max(identities):.3f}",
            "members": ",".join(members),
        })

    with PAIRWISE_TSV.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "cluster_id",
                "seq1",
                "seq2",
                "matches",
                "aligned_length",
                "identity_percent",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(all_rows)

    with SUMMARY_TSV.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "cluster_id",
                "member_count",
                "pair_count",
                "min_identity_percent",
                "mean_identity_percent",
                "max_identity_percent",
                "members",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print("Pairwise global-identity audit completed.")
    print(f"Wrote: {PAIRWISE_TSV}")
    print(f"Wrote: {SUMMARY_TSV}")
    print()
    print("Cluster summary:")
    for r in summary_rows:
        print(
            f"{r['cluster_id']}\t"
            f"n={r['member_count']}\t"
            f"pairs={r['pair_count']}\t"
            f"min={r['min_identity_percent']}\t"
            f"mean={r['mean_identity_percent']}\t"
            f"max={r['max_identity_percent']}"
        )


if __name__ == "__main__":
    main()
