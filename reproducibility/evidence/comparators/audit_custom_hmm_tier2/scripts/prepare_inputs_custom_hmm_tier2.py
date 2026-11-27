#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd


OUTDIR = Path("audit_custom_hmm_tier2")
INPUTS = OUTDIR / "inputs"
MANIFESTS = OUTDIR / "manifests"

CALIB_FEATURES = Path("./calibration_42_feature_table.csv")
REFSEQ_44_FASTA = Path("./halobacteria_44_HalC8_domains_extracted_by_hmmer.faa")
REFSEQ_18_SCORES = Path("./halobacteria_18_pipeline_scores.csv")

CALIB_OUT = INPUTS / "calibration_42_extracted_domains_from_feature_table.faa"
REFSEQ18_OUT = INPUTS / "refseq_18_calibration_set_novel_domains_from_44.faa"
REFSEQ18_MANIFEST = MANIFESTS / "refseq18_extraction_manifest.tsv"


def canon(x):
    s = str(x).strip().lstrip(">")
    m = re.search(r"(WP_\d+\.\d+)", s)
    if m:
        return m.group(1)
    return s.split()[0]


def write_wrapped(out, header, seq):
    out.write(f">{header}\n")
    for i in range(0, len(seq), 60):
        out.write(seq[i:i+60] + "\n")


def read_fasta(path):
    records = {}
    header = None
    seq = []

    def flush():
        if header is None:
            return
        acc = canon(header)
        records[acc] = (header.lstrip(">"), "".join(seq).replace("-", ""))

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                flush()
                header = line
                seq = []
            else:
                seq.append(line)
        flush()

    return records


def infer_id_column(df):
    for col in df.columns:
        if df[col].astype(str).str.contains(r"WP_\d+\.\d+", regex=True).any():
            return col
    for col in ["accession", "protein_accession", "seq_id", "id"]:
        if col in df.columns:
            return col
    raise SystemExit(f"Could not infer accession column. Columns: {list(df.columns)}")


def main():
    INPUTS.mkdir(parents=True, exist_ok=True)
    MANIFESTS.mkdir(parents=True, exist_ok=True)

    # Calibration n=42 domain FASTA from mature_domain_seq
    df = pd.read_csv(CALIB_FEATURES)
    required = {"accession", "mature_domain_seq"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns in {CALIB_FEATURES}: {missing}")

    with open(CALIB_OUT, "w") as out:
        for _, row in df.iterrows():
            acc = str(row["accession"]).strip()
            seq = str(row["mature_domain_seq"]).strip().replace("-", "")
            if not seq or seq.lower() == "nan":
                raise SystemExit(f"Empty mature_domain_seq for {acc}")
            write_wrapped(out, acc, seq)

    # RefSeq n=18 FASTA from 44 HMM-extracted domains
    records = read_fasta(REFSEQ_44_FASTA)
    scores = pd.read_csv(REFSEQ_18_SCORES)
    id_col = infer_id_column(scores)
    wanted = list(dict.fromkeys(canon(x) for x in scores[id_col].tolist()))

    missing = []
    with open(REFSEQ18_OUT, "w") as out:
        for acc in wanted:
            if acc not in records:
                missing.append(acc)
                continue
            header, seq = records[acc]
            write_wrapped(out, header, seq)

    with open(REFSEQ18_MANIFEST, "w") as m:
        m.write("accession\tstatus\n")
        for acc in wanted:
            m.write(f"{acc}\t{'found' if acc in records else 'missing'}\n")

    print("calibration_written:", CALIB_OUT)
    print("calibration_n:", sum(1 for line in open(CALIB_OUT) if line.startswith(">")))
    print("refseq18_written:", REFSEQ18_OUT)
    print("refseq18_n:", sum(1 for line in open(REFSEQ18_OUT) if line.startswith(">")))
    print("refseq18_missing:", ",".join(missing) if missing else "none")


if __name__ == "__main__":
    main()

