#!/usr/bin/env python3
"""Post hoc reconstruction, NOT the original v2.2 extractor.
Scope: the frozen single-domain-per-target RefSeq deployment.
Input: NR100 protein FASTA and retained HMMER hmmsearch domtblout.
Output: FASTA with accession and 1-based inclusive extraction coordinates.
"""
import argparse
from pathlib import Path
from math import isfinite


def fasta_records(handle):
    name = None
    chunks = []
    for line in handle:
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def extract(proteins, domains):
    hits = {}
    for line in domains:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split(maxsplit=22)
        acc = fields[0]
        n, start, end = map(int, (fields[2], fields[19], fields[20]))
        e = float(fields[12])
        if not (isfinite(e) and 0 <= e <= 1e-10):
            raise ValueError(f"Unretained domain: {acc}")
        if acc in hits or not (1 <= start <= end <= n):
            raise ValueError(f"Duplicate target or invalid envelope: {acc}")
        hits[acc] = (n, max(1, start - 2) if n > 120 else start, end)
    if not hits:
        raise ValueError("No retained domains")
    found = {}
    for acc, seq in fasta_records(proteins):
        if acc not in hits:
            continue
        n, start, end = hits[acc]
        if acc in found or len(seq) != n:
            raise ValueError(f"Duplicate target or length mismatch: {acc}")
        found[acc] = seq[start-1:end]
    if set(found) != set(hits):
        raise ValueError(f"Missing targets: {sorted(set(hits)-set(found))}")
    return {f"{acc}_domain_{start}_{end}": found[acc]
            for acc, (_, start, end) in hits.items()}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--proteins", type=Path, required=True)
    p.add_argument("--domtblout", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    with args.proteins.open() as proteins, args.domtblout.open() as domains:
        records = extract(proteins, domains)
    with args.output.open("x") as out:
        for name, seq in records.items():
            out.write(f">{name}\n{seq}\n")
    print(f"Extracted {len(records)} domains")


if __name__ == "__main__":
    main()
