#!/usr/bin/env bash
set -euo pipefail

source audit_custom_hmm_tier2/paths.env

mkdir -p audit_custom_hmm_tier2/{hmmer,logs,manifests}

shasum -a 256 "$HMM" "$CALIB_FASTA" "$CALIB_LABELS" "$CALIB_TIER2" "$REFSEQ_FASTA" "$REFSEQ_TIER2" \
  > audit_custom_hmm_tier2/manifests/input_sha256.txt

{
  echo -e "dataset\tfasta\tsequence_count"
  echo -e "calibration_42\t$CALIB_FASTA\t$(grep -c '^>' "$CALIB_FASTA")"
  echo -e "refseq_18\t$REFSEQ_FASTA\t$(grep -c '^>' "$REFSEQ_FASTA")"
} > audit_custom_hmm_tier2/manifests/input_sequence_counts.tsv

cat audit_custom_hmm_tier2/manifests/input_sequence_counts.tsv

hmmsearch \
  --max \
  -E 1000 \
  --domE 1000 \
  --domtblout audit_custom_hmm_tier2/hmmer/calibration_42_vs_HalC8_custom_v2.2.domtblout \
  "$HMM" \
  "$CALIB_FASTA" \
  > audit_custom_hmm_tier2/logs/calibration_42_vs_HalC8_custom_v2.2.hmmsearch.log

hmmsearch \
  --max \
  -E 1000 \
  --domE 1000 \
  --domtblout audit_custom_hmm_tier2/hmmer/refseq_18_vs_HalC8_custom_v2.2.domtblout \
  "$HMM" \
  "$REFSEQ_FASTA" \
  > audit_custom_hmm_tier2/logs/refseq_18_vs_HalC8_custom_v2.2.hmmsearch.log

ls -lh audit_custom_hmm_tier2/hmmer/*.domtblout
