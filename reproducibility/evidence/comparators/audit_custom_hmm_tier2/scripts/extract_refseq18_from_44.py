#!/usr/bin/env python3

import re
import sys
from pathlib import Path
import pandas as pd


SOURCE_FASTA = Path("./halobacteria_44_HalC8_domains_extracted_by_hmmer.faa")
SCORES_CSV = Path("./halobacteria_18_pipeline_scores.csv")
OUT_FASTA = Path("audit_custom_hmm_tier2/inputs/refseq_18_calibration_set_novel_domains_from_44.faa")
OUT_MANIFEST = Path("audit_custom_hmm_tier2/manifests/refseq18_extraction_manifest.tsv")


def canon(s):
    s = str(s).strip()
    m = re.search(r"(WP_\d+\.\d+)", s)
    if m:
        return m.group(1)
    return s.split()[0].lstrip(">")


def read_fasta(path):
    records = {}
    current_id = None
    header = None
    seq = []

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = (header, "".join(seq))
                header = line
                current_id = canon(line)
                seq = []
            else:
                seq.append(line.strip())

    if current_id is not None:
        records[current_id] = (header, "".join(seq))

    return records


def infer_id_column(df):
    candidates = [
        "accession",
        "protein_accession",
        "seq_id",
        "sequence_id",
        "id",
        "target",
        "name",
        "domain_id",
        "query",
    ]
    normalized = {re.sub(r"[^a-z0-9]+", "_", c.lower()).strip("_"): c for c in df.columns}

    for c in candidates:
        key = re.sub(r"[^a-z0-9]+", "_", c.lower()).strip("_")
        if key in normalized:
            return normalized[key]

    # Fallback: first column containing WP_ in at least one row.
    for col in df.columns:
        if df[col].astype(str).str.contains(r"WP_\d+\.\d+", regex=True).any():
            return col

    raise SystemExit(f"Could not infer ID column. Columns: {list(df.columns)}")


def main():
    if not SOURCE_FASTA.exists():
        raise SystemExit(f"Missing source FASTA: {SOURCE_FASTA}")
    if not SCORES_CSV.exists():
        raise SystemExit(f"Missing scores CSV: {SCORES_CSV}")

    records = read_fasta(SOURCE_FASTA)
    df = pd.read_csv(SCORES_CSV)
    id_col = infer_id_column(df)

    wanted = [canon(x) for x in df[id_col].tolist()]
    wanted_unique = list(dict.fromkeys(wanted))

    OUT_FASTA.parent.mkdir(parents=True, exist_ok=True)
    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)

    missing = []
    with open(OUT_FASTA, "w") as out:
        for acc in wanted_unique:
            if acc not in records:
                missing.append(acc)
                continue
            header, seq = records[acc]
            out.write(f"{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

    with open(OUT_MANIFEST, "w") as m:
        m.write("accession\tstatus\n")
        for acc in wanted_unique:
            m.write(f"{acc}\t{'found' if acc in records else 'missing'}\n")

    print("score_rows:", len(df))
    print("unique_score_ids:", len(wanted_unique))
    print("source_fasta_records:", len(records))
    print("written_records:", len(wanted_unique) - len(missing))
    print("missing:", ",".join(missing) if missing else "none")
    print("output:", OUT_FASTA)


if __name__ == "__main__":
    main()
