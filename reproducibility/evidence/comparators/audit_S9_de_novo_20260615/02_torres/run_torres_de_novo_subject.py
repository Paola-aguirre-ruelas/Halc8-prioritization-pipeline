from pathlib import Path
import csv, json, hashlib, subprocess, sys

ROOT = Path("audit_S9_de_novo_20260615")
OUT = ROOT / "02_torres"
OUT.mkdir(parents=True, exist_ok=True)

QUERY = Path("torres2025_archaeasins_DataS2.fasta")
TARGET_ORIGINAL = Path("audit_S8_de_novo_20260615/01_build_sets/calibration_39_CANONICAL_for_HMMER.faa")

def fail(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

def sha256_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def sha256_seq(s):
    return hashlib.sha256(s.encode()).hexdigest()

def read_fasta(path):
    records = []
    header = None
    seq = []
    for line in path.read_text(errors="ignore").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq).replace(" ", "").replace("\t", "").upper()))
            header = line[1:]
            seq = []
        else:
            seq.append(line)
    if header is not None:
        records.append((header, "".join(seq).replace(" ", "").replace("\t", "").upper()))
    return records

def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        fail("Command failed:\n" + " ".join(cmd) + "\nSTDERR:\n" + p.stderr)
    return p

for p in [QUERY, TARGET_ORIGINAL]:
    if not p.exists():
        fail(f"Missing required file: {p}")

queries = read_fasta(QUERY)
targets = read_fasta(TARGET_ORIGINAL)

if len(queries) != 265:
    fail(f"Expected 265 Torres DataS2 queries, observed {len(queries)}")
if len(targets) != 39:
    fail(f"Expected 39 HalC8 target sequences, observed {len(targets)}")

query_lengths = [len(s) for _, s in queries]
target_lengths = [len(s) for _, s in targets]

target_sanitized = OUT / "HalC8_n39_CANONICAL_for_BLAST_subject_sanitized.faa"
target_map = OUT / "HalC8_n39_CANONICAL_for_BLAST_subject_sanitized_id_map.tsv"

with target_sanitized.open("w") as f_faa, target_map.open("w", newline="") as f_map:
    w = csv.DictWriter(
        f_map,
        fieldnames=["subject_id", "original_id", "original_header", "length", "sequence_sha256"],
        delimiter="\t",
    )
    w.writeheader()
    for i, (header, seq) in enumerate(targets, start=1):
        original_id = header.split("|")[0].strip().split()[0]
        subject_id = f"HalC8_{i:02d}"
        f_faa.write(f">{subject_id}\n{seq}\n")
        w.writerow({
            "subject_id": subject_id,
            "original_id": original_id,
            "original_header": header,
            "length": len(seq),
            "sequence_sha256": sha256_seq(seq),
        })

query_audit = OUT / "torres_DataS2_query_audit.tsv"
with query_audit.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["query_id", "length", "sequence_sha256", "header"], delimiter="\t")
    w.writeheader()
    for h, s in queries:
        w.writerow({
            "query_id": h.split()[0],
            "length": len(s),
            "sequence_sha256": sha256_seq(s),
            "header": h,
        })

outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"

strict_out = OUT / "torres_vs_HalC8_n39_blastp_strict_subject.tsv"
short_out = OUT / "torres_vs_HalC8_n39_blastp_short_subject.tsv"

strict_cmd = [
    "blastp",
    "-query", str(QUERY),
    "-subject", str(target_sanitized),
    "-out", str(strict_out),
    "-outfmt", outfmt,
    "-evalue", "1e-3",
]

short_cmd = [
    "blastp",
    "-task", "blastp-short",
    "-query", str(QUERY),
    "-subject", str(target_sanitized),
    "-out", str(short_out),
    "-outfmt", outfmt,
    "-evalue", "10",
]

run(strict_cmd)
run(short_cmd)

FIELDS = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen","qcovs"]

def read_blast(path):
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with path.open() as f:
        for line in f:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            r = dict(zip(FIELDS, vals))
            r["pident"] = float(r["pident"])
            r["length"] = int(r["length"])
            r["evalue"] = float(r["evalue"])
            r["bitscore"] = float(r["bitscore"])
            r["qlen"] = int(r["qlen"])
            r["slen"] = int(r["slen"])
            r["qcovs"] = float(r["qcovs"])
            rows.append(r)
    return rows

strict_rows = read_blast(strict_out)
short_rows = read_blast(short_out)

strict_hits = [r for r in strict_rows if r["pident"] >= 25.0]

short_meaningful = [
    r for r in short_rows
    if r["pident"] >= 25.0 and r["qcovs"] >= 50.0 and r["length"] >= 8
]

best_short_by_query = {}
for r in short_meaningful:
    q = r["qseqid"]
    old = best_short_by_query.get(q)
    if old is None or (r["evalue"], -r["bitscore"], -r["pident"]) < (old["evalue"], -old["bitscore"], -old["pident"]):
        best_short_by_query[q] = r

short_meaningful_best = list(best_short_by_query.values())

strict_hits_path = OUT / "torres_vs_HalC8_n39_strict_subject_hits_filtered.tsv"
short_hits_path = OUT / "torres_vs_HalC8_n39_short_subject_meaningful_hits_filtered.tsv"
summary_tsv = OUT / "torres_vs_HalC8_n39_CANONICAL_subject_summary.tsv"
manifest_json = OUT / "torres_vs_HalC8_n39_CANONICAL_subject_manifest.json"

def write_rows(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in FIELDS})

write_rows(strict_hits_path, strict_hits)
write_rows(short_hits_path, short_meaningful_best)

summary = {
    "analysis": "Torres2025 archaeasin DataS2 vs canonical HalC8 n39 mature-domain set using blastp -subject",
    "n_queries": len(queries),
    "query_length_min": min(query_lengths),
    "query_length_max": max(query_lengths),
    "query_length_mean": round(sum(query_lengths) / len(query_lengths), 3),
    "n_targets": len(targets),
    "target_length_min": min(target_lengths),
    "target_length_max": max(target_lengths),
    "strict_blastp_raw_rows": len(strict_rows),
    "strict_blastp_hits_E_le_1e_3_identity_ge_25": len(strict_hits),
    "short_blastp_raw_rows": len(short_rows),
    "short_blastp_meaningful_best_hits_identity_ge_25_qcov_ge_50_aln_ge_8": len(short_meaningful_best),
}

with summary_tsv.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary.keys()), delimiter="\t")
    w.writeheader()
    w.writerow(summary)

manifest = {
    "analysis": "S9 Part B de novo Torres2025 archaeasin complementarity audit; BLASTP subject mode used because makeblastdb failed locally with memory-map errors",
    "query_fasta": str(QUERY),
    "query_fasta_sha256": sha256_file(QUERY),
    "target_fasta_original": str(TARGET_ORIGINAL),
    "target_fasta_original_sha256": sha256_file(TARGET_ORIGINAL),
    "target_fasta_subject_sanitized": str(target_sanitized),
    "target_fasta_subject_sanitized_sha256": sha256_file(target_sanitized),
    "target_id_map": str(target_map),
    "target_id_map_sha256": sha256_file(target_map),
    "strict_blastp_command": " ".join(strict_cmd),
    "short_blastp_command": " ".join(short_cmd),
    "strict_hit_definition": "blastp -subject, E<=1e-3, identity>=25%",
    "short_meaningful_hit_definition": "blastp-short -subject, E<=10 reported hits post-filtered by identity>=25%, query coverage>=50%, alignment length>=8 aa; best hit per query retained",
    "outputs": {
        "query_audit": str(query_audit),
        "target_sanitized": str(target_sanitized),
        "target_map": str(target_map),
        "strict_raw": str(strict_out),
        "short_raw": str(short_out),
        "strict_filtered": str(strict_hits_path),
        "short_filtered": str(short_hits_path),
        "summary": str(summary_tsv),
    },
    "sha256": {
        "query_audit": sha256_file(query_audit),
        "target_sanitized": sha256_file(target_sanitized),
        "target_map": sha256_file(target_map),
        "strict_raw": sha256_file(strict_out),
        "short_raw": sha256_file(short_out),
        "strict_filtered": sha256_file(strict_hits_path),
        "short_filtered": sha256_file(short_hits_path),
        "summary": sha256_file(summary_tsv),
    }
}

manifest_json.write_text(json.dumps(manifest, indent=2))

print("== Torres de novo subject-mode summary ==")
for k, v in summary.items():
    print(f"{k}: {v}")

print("\n== Filtered strict hits ==")
print(strict_hits_path.read_text())

print("== Filtered short-peptide meaningful hits ==")
print(short_hits_path.read_text())

print("== WROTE ==")
for p in [target_sanitized, target_map, query_audit, strict_out, short_out, strict_hits_path, short_hits_path, summary_tsv, manifest_json]:
    print(p, sha256_file(p))
