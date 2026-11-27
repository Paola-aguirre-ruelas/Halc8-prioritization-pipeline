#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$BASE"

INPUT="input/halobacteria_13_score_positive.fasta"

mkdir -p cdhit mmseqs logs tables

echo "Canonical n=13 HalC8 score-positive clustering"
echo "Base directory: $BASE"
echo "Input FASTA: $INPUT"

if [[ ! -f "$INPUT" ]]; then
  echo "ERROR: input FASTA not found: $INPUT" >&2
  exit 1
fi

NSEQ="$(grep -c '^>' "$INPUT")"
echo "Input sequence count: $NSEQ"

if [[ "$NSEQ" != "13" ]]; then
  echo "ERROR: expected 13 sequences, found $NSEQ" >&2
  exit 1
fi

{
  echo "## Date"
  date
  echo
  echo "## Base directory"
  pwd
  echo
  echo "## Input checksum"
  shasum -a 256 "$INPUT"
  echo
  echo "## FASTA headers"
  grep '^>' "$INPUT"
  echo
  echo "## Tool paths"
  command -v cd-hit
  command -v mmseqs
  command -v needleall || true
  echo
  echo "## Tool versions"
  cd-hit -h 2>&1 | head -5 || true
  mmseqs version 2>&1 || true
  needleall -help 2>&1 | head -5 || true
} > logs/run_environment.txt

for THR in 1.00 0.95 0.90 0.80; do
  echo "=== CD-HIT threshold $THR ==="
  cd-hit \
    -i "$INPUT" \
    -o "cdhit/cdhit_n13_c${THR}.fasta" \
    -c "$THR" \
    -n 5 \
    -G 1 \
    -aS 0.8 \
    -g 1 \
    -d 0 \
    -T 4 \
    > "logs/cdhit_n13_c${THR}.log" 2>&1

  echo "=== MMseqs2 low-level threshold $THR ==="

  MMDB="mmseqs/mmseqs_n13_c${THR}_db"
  CLU="mmseqs/mmseqs_n13_c${THR}_clu"
  TMP="mmseqs/tmp_c${THR}"
  TSV="mmseqs/mmseqs_n13_c${THR}_cluster.tsv"

  rm -rf "$TMP" "$MMDB"* "$CLU"* "$TSV"

  {
    echo "## mmseqs createdb"
    mmseqs createdb "$INPUT" "$MMDB"

    echo
    echo "## mmseqs cluster"
    mmseqs cluster \
      "$MMDB" \
      "$CLU" \
      "$TMP" \
      --min-seq-id "$THR" \
      -c 0.8 \
      --cov-mode 0 \
      --cluster-mode 0 \
      -s 7.5 \
      --threads 4

    echo
    echo "## mmseqs createtsv"
    mmseqs createtsv "$MMDB" "$MMDB" "$CLU" "$TSV"
  } > "logs/mmseqs_n13_c${THR}.log" 2>&1
done

echo "DONE: primary CD-HIT/MMseqs2 clustering completed."
