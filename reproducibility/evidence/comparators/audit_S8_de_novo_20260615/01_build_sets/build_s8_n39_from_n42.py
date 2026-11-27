from pathlib import Path
import csv, json, hashlib, subprocess, sys
from collections import Counter

ROOT = Path("audit_S8_de_novo_20260615")
INVENTORY = ROOT / "00_input_inventory"
BUILD = ROOT / "01_build_sets"
LOGS = ROOT / "logs"

BUILD.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)

N42_FASTA = Path("calibration_42_for_amp_predictors.fasta")
LABELS_CSV = Path("calibration_42_labels.csv")
CUSTOM_HMM = Path("HalC8_custom_v2.2.hmm")

EXCLUDED_FROM_N39 = {
    "AGB32219.1": "fragment/domain-unstable calibration sequence; excluded from mature-domain HMM benchmark",
    "WP_138653937.1": "fragment/domain-unstable calibration sequence; excluded from mature-domain HMM benchmark",
    "WP_211249102.1": "fragment/domain-unstable calibration sequence; excluded from mature-domain HMM benchmark",
}

DOMAIN_COMPATIBLE_MAX_LEN = 130

def fail(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

def sha256_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def sha256_seq(seq):
    return hashlib.sha256(seq.encode()).hexdigest()

def norm_acc(header):
    h = header[1:] if header.startswith(">") else header
    return h.split("|")[0].strip().split()[0]

def read_fasta(path):
    records = {}
    order = []
    header = None
    seq = []
    for line in path.read_text(errors="ignore").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                acc = norm_acc(header)
                sequence = "".join(seq).replace(" ", "").replace("\t", "")
                if acc in records:
                    fail(f"duplicate FASTA accession: {acc}")
                records[acc] = {"header": header[1:], "sequence": sequence}
                order.append(acc)
            header = line
            seq = []
        else:
            seq.append(line)
    if header is not None:
        acc = norm_acc(header)
        sequence = "".join(seq).replace(" ", "").replace("\t", "")
        if acc in records:
            fail(f"duplicate FASTA accession: {acc}")
        records[acc] = {"header": header[1:], "sequence": sequence}
        order.append(acc)
    return records, order

def read_labels(path):
    rows = list(csv.DictReader(path.open()))
    required = {"accession", "binary_label", "tfm_verdict"}
    if not required.issubset(rows[0].keys()):
        fail(f"labels file must contain {required}; observed {rows[0].keys()}")

    labels = {}
    for r in rows:
        acc = r["accession"].strip()
        labels[acc] = {
            "binary_label": r["binary_label"].strip(),
            "tfm_verdict": r["tfm_verdict"].strip(),
        }
    return labels

def run_hmmer_for_domain_extraction():
    tbl = BUILD / "n42_vs_HalC8_custom_for_domain_extraction.tblout"
    domtbl = BUILD / "n42_vs_HalC8_custom_for_domain_extraction.domtblout"
    stdout = BUILD / "n42_vs_HalC8_custom_for_domain_extraction.hmmsearch.txt"

    cmd = [
        "hmmsearch",
        "--cpu", "4",
        "--noali",
        "--tblout", str(tbl),
        "--domtblout", str(domtbl),
        str(CUSTOM_HMM),
        str(N42_FASTA),
    ]

    with stdout.open("w") as out:
        p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)

    if p.returncode != 0:
        fail(f"hmmsearch failed:\n{p.stderr}")

    return cmd, tbl, domtbl, stdout

def parse_best_domains(domtbl):
    best = {}
    with domtbl.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            target = parts[0]
            full_e = float(parts[6])
            full_score = float(parts[7])
            c_e = float(parts[11])
            i_e = float(parts[12])
            dom_score = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            env_from = int(parts[19])
            env_to = int(parts[20])
            acc = float(parts[21])

            hit = {
                "full_evalue": full_e,
                "full_score": full_score,
                "conditional_evalue": c_e,
                "independent_evalue": i_e,
                "domain_score": dom_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "env_from": env_from,
                "env_to": env_to,
                "alignment_accuracy": acc,
            }

            if target not in best:
                best[target] = hit
            else:
                old = best[target]
                if (hit["independent_evalue"], -hit["domain_score"]) < (old["independent_evalue"], -old["domain_score"]):
                    best[target] = hit
    return best

for p in [N42_FASTA, LABELS_CSV, CUSTOM_HMM]:
    if not p.exists():
        fail(f"missing required input: {p}")

records, order = read_fasta(N42_FASTA)
labels = read_labels(LABELS_CSV)

print("== Input audit ==")
print("n42_fasta:", N42_FASTA, "records:", len(records), "sha256:", sha256_file(N42_FASTA))
print("labels:", LABELS_CSV, "rows:", len(labels), "sha256:", sha256_file(LABELS_CSV))
print("custom_hmm:", CUSTOM_HMM, "sha256:", sha256_file(CUSTOM_HMM))

missing_labels = sorted(set(records) - set(labels))
missing_sequences = sorted(set(labels) - set(records))
if missing_labels:
    fail(f"FASTA records missing labels: {missing_labels}")
if missing_sequences:
    fail(f"label rows missing FASTA sequence: {missing_sequences}")
if len(records) != 42 or len(labels) != 42:
    fail(f"expected 42 FASTA records and 42 labels; observed {len(records)} FASTA / {len(labels)} labels")

print("label_counts_n42:", dict(Counter(v["binary_label"] for v in labels.values())))

hmmer_cmd, tbl, domtbl, hmmer_stdout = run_hmmer_for_domain_extraction()
best_domains = parse_best_domains(domtbl)

audit42_path = BUILD / "calibration_42_CANONICAL_input_audit.tsv"
eligibility_path = BUILD / "calibration_42_to_39_CANONICAL_eligibility.tsv"
n39_fasta = BUILD / "calibration_39_CANONICAL_for_HMMER.faa"
n39_labels = BUILD / "calibration_39_CANONICAL_labels.tsv"
manifest_path = BUILD / "calibration_39_CANONICAL_build_manifest.json"

audit_rows = []
elig_rows = []
n39_rows = []

for acc in order:
    seq = records[acc]["sequence"]
    lab = labels[acc]
    length = len(seq)
    cys = seq.count("C")
    seq_hash = sha256_seq(seq)

    audit_rows.append({
        "accession": acc,
        "binary_label": lab["binary_label"],
        "tfm_verdict": lab["tfm_verdict"],
        "original_length": length,
        "original_cys": cys,
        "original_sequence_sha256": seq_hash,
        "original_header": records[acc]["header"],
    })

    if acc in EXCLUDED_FROM_N39:
        elig_rows.append({
            "accession": acc,
            "binary_label": lab["binary_label"],
            "tfm_verdict": lab["tfm_verdict"],
            "original_length": length,
            "included_in_n39": "NO",
            "domain_source": "excluded",
            "domain_start": "",
            "domain_end": "",
            "domain_length": "",
            "domain_cys": "",
            "domain_i_evalue": "",
            "exclusion_reason": EXCLUDED_FROM_N39[acc],
        })
        continue

    if length <= DOMAIN_COMPATIBLE_MAX_LEN:
        domain_seq = seq
        source = "full_sequence_domain_compatible"
        start = 1
        end = length
        i_e = ""
    else:
        if acc not in best_domains:
            fail(f"long sequence {acc} length {length} has no HalC8_custom HMM domain hit for extraction")
        hit = best_domains[acc]
        start = hit["ali_from"]
        end = hit["ali_to"]
        domain_seq = seq[start-1:end]
        source = "extracted_best_HalC8_custom_domain"
        i_e = hit["independent_evalue"]

    n39_rows.append({
        "accession": acc,
        "binary_label": lab["binary_label"],
        "tfm_verdict": lab["tfm_verdict"],
        "domain_source": source,
        "original_length": length,
        "domain_start": start,
        "domain_end": end,
        "domain_length": len(domain_seq),
        "domain_cys": domain_seq.count("C"),
        "domain_i_evalue": i_e,
        "domain_sequence_sha256": sha256_seq(domain_seq),
        "domain_sequence": domain_seq,
    })

    elig_rows.append({
        "accession": acc,
        "binary_label": lab["binary_label"],
        "tfm_verdict": lab["tfm_verdict"],
        "original_length": length,
        "included_in_n39": "YES",
        "domain_source": source,
        "domain_start": start,
        "domain_end": end,
        "domain_length": len(domain_seq),
        "domain_cys": domain_seq.count("C"),
        "domain_i_evalue": i_e,
        "exclusion_reason": "",
    })

counts39 = Counter(r["binary_label"] for r in n39_rows)

if len(n39_rows) != 39:
    fail(f"expected n39=39, observed {len(n39_rows)}")
if len({r["accession"] for r in n39_rows}) != 39:
    fail("duplicate accessions in n39")
if counts39.get("POSITIVE", 0) != 23 or counts39.get("NEGATIVE", 0) != 16:
    fail(f"expected n39 label counts POSITIVE=23 NEGATIVE=16; observed {dict(counts39)}")

with audit42_path.open("w", newline="") as f:
    fieldnames = list(audit_rows[0].keys())
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(audit_rows)

with eligibility_path.open("w", newline="") as f:
    fieldnames = list(elig_rows[0].keys())
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(elig_rows)

with n39_labels.open("w", newline="") as f:
    fieldnames = [
        "accession", "binary_label", "tfm_verdict", "domain_source",
        "original_length", "domain_start", "domain_end", "domain_length",
        "domain_cys", "domain_i_evalue", "domain_sequence_sha256"
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    for r in n39_rows:
        w.writerow({k: r[k] for k in fieldnames})

with n39_fasta.open("w") as f:
    for r in n39_rows:
        f.write(
            f">{r['accession']}|{r['binary_label']}|{r['tfm_verdict']}|"
            f"{r['domain_source']}|origlen{r['original_length']}|"
            f"dom{r['domain_start']}-{r['domain_end']}|len{r['domain_length']}|"
            f"cys{r['domain_cys']}|seqsha256:{r['domain_sequence_sha256']}\n"
        )
        f.write(r["domain_sequence"] + "\n")

manifest = {
    "analysis": "S8 de novo canonical rebuild of n=39 mature-domain-compatible HMM benchmark set from n=42 calibration universe",
    "n42_fasta": str(N42_FASTA),
    "n42_fasta_sha256": sha256_file(N42_FASTA),
    "labels_csv": str(LABELS_CSV),
    "labels_csv_sha256": sha256_file(LABELS_CSV),
    "custom_hmm_for_domain_extraction": str(CUSTOM_HMM),
    "custom_hmm_sha256": sha256_file(CUSTOM_HMM),
    "hmmer_domain_extraction_command": " ".join(hmmer_cmd),
    "n42_records": len(records),
    "n42_label_counts": dict(Counter(v["binary_label"] for v in labels.values())),
    "excluded_from_n39": EXCLUDED_FROM_N39,
    "n39_records": len(n39_rows),
    "n39_label_counts": dict(counts39),
    "n39_domain_source_counts": dict(Counter(r["domain_source"] for r in n39_rows)),
    "n39_domain_length_range": [min(r["domain_length"] for r in n39_rows), max(r["domain_length"] for r in n39_rows)],
    "outputs": {
        "audit42": str(audit42_path),
        "eligibility": str(eligibility_path),
        "n39_fasta": str(n39_fasta),
        "n39_labels": str(n39_labels),
        "domain_extraction_tblout": str(tbl),
        "domain_extraction_domtblout": str(domtbl),
        "domain_extraction_hmmsearch_txt": str(hmmer_stdout),
    },
}
manifest_path.write_text(json.dumps(manifest, indent=2))

print("\n== n39 canonical build summary ==")
print("n39_records:", len(n39_rows))
print("n39_label_counts:", dict(counts39))
print("n39_domain_source_counts:", dict(Counter(r["domain_source"] for r in n39_rows)))
print("n39_domain_length_range:", manifest["n39_domain_length_range"])
print("excluded_from_n39:", sorted(EXCLUDED_FROM_N39))

print("\n== WROTE ==")
for p in [audit42_path, eligibility_path, n39_fasta, n39_labels, manifest_path, tbl, domtbl, hmmer_stdout]:
    print(p, sha256_file(p))
