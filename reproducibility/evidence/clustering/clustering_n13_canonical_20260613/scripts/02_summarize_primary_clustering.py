#!/usr/bin/env python3

from pathlib import Path
import csv
import re
from collections import defaultdict

SCRIPT_DIR = Path(__file__).resolve().parent
BASE = SCRIPT_DIR.parent
THRESHOLDS = ["1.00", "0.95", "0.90", "0.80"]

def wp_id(text: str) -> str:
    m = re.search(r"(WP_\d+\.\d+)", text)
    return m.group(1) if m else text.strip()

def parse_cdhit_clstr(path: Path):
    clusters = []
    current = []

    with path.open() as f:
        for raw in f:
            line = raw.strip()
            if line.startswith(">Cluster"):
                if current:
                    clusters.append(current)
                current = []
                continue

            if ">" in line:
                token = line.split(">", 1)[1].split("...", 1)[0]
                current.append(wp_id(token))

    if current:
        clusters.append(current)

    return clusters

def parse_mmseqs_cluster_tsv(path: Path):
    clusters_by_rep = defaultdict(list)

    with path.open() as f:
        for raw in f:
            if not raw.strip():
                continue
            rep, member = raw.rstrip("\n").split("\t")[:2]
            clusters_by_rep[wp_id(rep)].append(wp_id(member))

    return list(clusters_by_rep.values())

def norm_cluster_set(clusters):
    return sorted([tuple(sorted(c)) for c in clusters])

rows = []

for thr in THRESHOLDS:
    cdhit_file = BASE / "cdhit" / f"cdhit_n13_c{thr}.fasta.clstr"
    mmseqs_file = BASE / "mmseqs" / f"mmseqs_n13_c{thr}_cluster.tsv"

    if not cdhit_file.exists():
        raise FileNotFoundError(f"Missing CD-HIT file: {cdhit_file}")

    if not mmseqs_file.exists():
        raise FileNotFoundError(f"Missing MMseqs2 cluster TSV: {mmseqs_file}")

    cdhit_clusters = parse_cdhit_clstr(cdhit_file)
    mmseqs_clusters = parse_mmseqs_cluster_tsv(mmseqs_file)

    rows.append({
        "threshold": thr,
        "cdhit_clusters": len(cdhit_clusters),
        "mmseqs_clusters": len(mmseqs_clusters),
        "concordant_count": len(cdhit_clusters) == len(mmseqs_clusters),
        "identical_membership": norm_cluster_set(cdhit_clusters) == norm_cluster_set(mmseqs_clusters),
    })

tables_dir = BASE / "tables"
tables_dir.mkdir(exist_ok=True)

counts_csv = tables_dir / "cluster_counts_n13_primary.csv"
with counts_csv.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "threshold",
            "cdhit_clusters",
            "mmseqs_clusters",
            "concordant_count",
            "identical_membership",
        ],
    )
    writer.writeheader()
    writer.writerows(rows)

clusters90 = parse_cdhit_clstr(BASE / "cdhit" / "cdhit_n13_c0.90.fasta.clstr")

membership_tsv = tables_dir / "cluster_membership_n13_cdhit_0.90.tsv"
with membership_tsv.open("w") as f:
    f.write("cluster_id\tmember_count\tmembers\n")
    for i, members in enumerate(clusters90, start=1):
        f.write(f"C{i}\t{len(members)}\t{','.join(members)}\n")

print("Cluster counts:")
for r in rows:
    print(
        f"{r['threshold']}\t"
        f"CD-HIT={r['cdhit_clusters']}\t"
        f"MMseqs2={r['mmseqs_clusters']}\t"
        f"count_concordant={r['concordant_count']}\t"
        f"membership_identical={r['identical_membership']}"
    )

print("\nCD-HIT 0.90 membership:")
for i, members in enumerate(clusters90, start=1):
    print(f"C{i}\tn={len(members)}\t{', '.join(members)}")

print(f"\nWrote: {counts_csv}")
print(f"Wrote: {membership_tsv}")
