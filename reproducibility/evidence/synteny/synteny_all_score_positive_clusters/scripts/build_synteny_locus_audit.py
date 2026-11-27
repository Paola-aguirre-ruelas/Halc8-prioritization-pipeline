#!/usr/bin/env python3

from pathlib import Path
import csv
import re
from collections import defaultdict

ALIGNMENTS = Path("results/reference_vs_candidates_alignments.tsv")
WINDOW_MANIFEST = Path("metadata/window_extraction_manifest_with_eval_status.tsv")
OUT = Path("results/halc8_synteny_locus_audit.tsv")

REF_GENE_MAP = {
    "AAZ66750.1": "halU",
    "AAT67988.1": "halR",
    "AAQ82548.1": "proC8",
    "AAT67989.2": "halT1",
    "AAZ66751.1": "halT2",
    "AAZ66752.1": "halT3",
}

GENES = ["halU", "halR", "proC8", "halT1", "halT2", "halT3"]


def load_manifest(path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            clean = {k.strip(): (v.strip() if isinstance(v, str) else v) for k, v in r.items()}
            rows.append(clean)
    return rows


def parse_reference_alignments(path):
    links = defaultdict(dict)
    current_locus = None

    with open(path) as f:
        for line in f:
            line = line.strip()

            if line.startswith("DQ137141.1_AS7092_HalC8_reference vs "):
                current_locus = line.split(" vs ", 1)[1].strip()
                continue

            if not current_locus:
                continue

            if not line or line.startswith("-") or line.startswith("Query"):
                continue

            parts = line.split("\t")
            if len(parts) != 4:
                continue

            query, target, identity, similarity = parts
            if query not in REF_GENE_MAP:
                continue

            gene = REF_GENE_MAP[query]
            links[current_locus][gene] = {
                "target": target,
                "identity": identity,
                "similarity": similarity,
            }

    return links


def locus_name_from_manifest_row(r):
    protein = r["protein_accession"].replace(".", "_")
    nuccore = r["nuccore_accession"].replace(".", "_")
    return f"{protein}__{nuccore}__pm15kb"


def classify(linked_genes, evaluable_status):
    linked = set(linked_genes)
    flanking = linked - {"proC8"}

    if evaluable_status == "LIMITED_EDGE_TRUNCATED":
        if "proC8" in linked and len(flanking) >= 3:
            return "limited_conserved", "edge-truncated window but proC8 plus >=3 flanking AS7092 genes detected"
        if len(linked) > 0:
            return "limited_partial", "edge-truncated window with some AS7092-linked genes detected"
        return "not_assessable", "edge-truncated window with no AS7092-linked genes detected"

    if "proC8" in linked and len(flanking) >= 3:
        return "conserved", "proC8 plus >=3 flanking AS7092 genes detected"

    if "proC8" in linked and 1 <= len(flanking) <= 2:
        return "partial", "proC8 plus 1-2 flanking AS7092 genes detected"

    if "proC8" not in linked and len(flanking) >= 3:
        return "partial", ">=3 flanking AS7092 genes detected without proC8 link"

    if len(linked) == 0:
        return "unsupported", "no AS7092-linked genes detected"

    return "weak_or_ambiguous", "limited AS7092-linked signal below predefined partial threshold"


def main():
    manifest_rows = load_manifest(WINDOW_MANIFEST)
    links = parse_reference_alignments(ALIGNMENTS)

    output_rows = []

    for r in manifest_rows:
        locus_name = locus_name_from_manifest_row(r)
        gene_links = links.get(locus_name, {})
        linked_genes = set(gene_links.keys())
        evaluable_status = r.get("synteny_evaluable_status", "")

        classification, reason = classify(linked_genes, evaluable_status)

        out = {
            "candidate_accession": r["protein_accession"],
            "nuccore_accession": r["nuccore_accession"],
            "evaluable_status": evaluable_status,
            "linked_halU": "YES" if "halU" in linked_genes else "NO",
            "linked_halR": "YES" if "halR" in linked_genes else "NO",
            "linked_proC8": "YES" if "proC8" in linked_genes else "NO",
            "linked_halT1": "YES" if "halT1" in linked_genes else "NO",
            "linked_halT2": "YES" if "halT2" in linked_genes else "NO",
            "linked_halT3": "YES" if "halT3" in linked_genes else "NO",
            "target_halU": gene_links.get("halU", {}).get("target", ""),
            "target_halR": gene_links.get("halR", {}).get("target", ""),
            "target_proC8": gene_links.get("proC8", {}).get("target", ""),
            "target_halT1": gene_links.get("halT1", {}).get("target", ""),
            "target_halT2": gene_links.get("halT2", {}).get("target", ""),
            "target_halT3": gene_links.get("halT3", {}).get("target", ""),
            "identity_halU": gene_links.get("halU", {}).get("identity", ""),
            "identity_halR": gene_links.get("halR", {}).get("identity", ""),
            "identity_proC8": gene_links.get("proC8", {}).get("identity", ""),
            "identity_halT1": gene_links.get("halT1", {}).get("identity", ""),
            "identity_halT2": gene_links.get("halT2", {}).get("identity", ""),
            "identity_halT3": gene_links.get("halT3", {}).get("identity", ""),
            "n_reference_genes_linked": len(linked_genes),
            "n_flanking_genes_linked_excluding_proC8": len(linked_genes - {"proC8"}),
            "classification": classification,
            "reason": reason,
        }

        output_rows.append(out)

    fieldnames = list(output_rows[0].keys())

    OUT.parent.mkdir(exist_ok=True)
    with open(OUT, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Wrote {OUT} with {len(output_rows)} rows")


if __name__ == "__main__":
    main()
