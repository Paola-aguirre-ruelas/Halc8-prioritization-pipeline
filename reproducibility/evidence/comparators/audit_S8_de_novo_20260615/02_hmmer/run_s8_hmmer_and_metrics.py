from pathlib import Path
import csv, json, hashlib, subprocess, sys, math
from collections import Counter

ROOT = Path("audit_S8_de_novo_20260615")
BUILD = ROOT / "01_build_sets"
HMMER = ROOT / "02_hmmer"
METRICS = ROOT / "03_metrics"
LOGS = ROOT / "logs"

HMMER.mkdir(parents=True, exist_ok=True)
METRICS.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)

FASTA = BUILD / "calibration_39_CANONICAL_for_HMMER.faa"
LABELS = BUILD / "calibration_39_CANONICAL_labels.tsv"

HMM_MODELS = {
    "HalC8_custom_v2.2": Path("HalC8_custom_v2.2.hmm"),
    "TIGR04449": Path("TIGR04449.hmm"),
    "PF22862": Path("PF22862.hmm"),
}

PRIMARY_IE_THRESHOLD = 1e-10
REPORT_E = "1000"
REPORT_DOME = "1000"

def fail(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

def sha256_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def get_hmmer_version():
    p = subprocess.run(["hmmsearch", "-h"], capture_output=True, text=True)
    if p.returncode != 0:
        fail("hmmsearch not available")
    for line in p.stdout.splitlines():
        if "HMMER" in line:
            return line.strip()
    return "unknown"

def read_labels(path):
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    labels = {}
    for r in rows:
        acc = r["accession"].strip()
        labels[acc] = {
            "binary_label": r["binary_label"].strip(),
            "tfm_verdict": r["tfm_verdict"].strip(),
            "length": int(r["domain_length"]),
            "cys": int(r["domain_cys"]),
            "domain_source": r["domain_source"],
        }
    return labels

def accession_from_target(target):
    return target.split("|")[0].strip()

def parse_tblout(path):
    # sequence-level hits
    hits = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            target = parts[0]
            acc = accession_from_target(target)
            evalue = float(parts[4])
            score = float(parts[5])
            bias = float(parts[6])
            if acc not in hits or evalue < hits[acc]["seq_evalue"]:
                hits[acc] = {
                    "target": target,
                    "seq_evalue": evalue,
                    "seq_score": score,
                    "seq_bias": bias,
                }
    return hits

def parse_domtblout(path):
    # best independent-domain hit per accession
    hits = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            target = parts[0]
            acc = accession_from_target(target)

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
            acc_prob = float(parts[21])

            hit = {
                "target": target,
                "full_evalue_domtbl": full_e,
                "full_score_domtbl": full_score,
                "conditional_evalue": c_e,
                "independent_evalue": i_e,
                "domain_score": dom_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "env_from": env_from,
                "env_to": env_to,
                "alignment_accuracy": acc_prob,
            }

            if acc not in hits:
                hits[acc] = hit
            else:
                old = hits[acc]
                if (hit["independent_evalue"], -hit["domain_score"]) < (old["independent_evalue"], -old["domain_score"]):
                    hits[acc] = hit
    return hits

def compute_metrics(labels, predicted_positive):
    tp = fp = tn = fn = 0
    for acc, info in labels.items():
        truth_pos = info["binary_label"] == "POSITIVE"
        pred_pos = acc in predicted_positive
        if truth_pos and pred_pos:
            tp += 1
        elif truth_pos and not pred_pos:
            fn += 1
        elif not truth_pos and pred_pos:
            fp += 1
        else:
            tn += 1

    n = tp + fp + tn + fn
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    f1 = 2 * precision * sens / (precision + sens) if precision + sens else float("nan")
    acc = (tp + tn) / n if n else float("nan")

    po = acc
    pe = (((tp + fp) / n) * ((tp + fn) / n) + ((tn + fn) / n) * ((tn + fp) / n)) if n else float("nan")
    kappa = (po - pe) / (1 - pe) if (1 - pe) else float("nan")

    denom = math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    mcc = ((tp*tn) - (fp*fn)) / denom if denom else float("nan")

    return {
        "n": n,
        "TP": tp, "FP": fp, "TN": tn, "FN": fn,
        "sensitivity": sens,
        "specificity": spec,
        "precision": precision,
        "F1": f1,
        "accuracy": acc,
        "kappa": kappa,
        "MCC": mcc,
    }

def fmt(x):
    if isinstance(x, float):
        if math.isnan(x):
            return "NA"
        return f"{x:.6g}"
    return x

for p in [FASTA, LABELS]:
    if not p.exists():
        fail(f"missing required input: {p}")

labels = read_labels(LABELS)
if len(labels) != 39:
    fail(f"expected 39 labels, observed {len(labels)}")

if Counter(v["binary_label"] for v in labels.values()) != Counter({"POSITIVE": 23, "NEGATIVE": 16}):
    fail(f"unexpected label counts: {Counter(v['binary_label'] for v in labels.values())}")

hmmer_version = get_hmmer_version()
print("== HMMER version ==")
print(hmmer_version)

all_per_sequence = []
summary_rows = []
commands = {}

for model_name, hmm_path in HMM_MODELS.items():
    if not hmm_path.exists():
        fail(f"missing HMM: {hmm_path}")

    prefix = HMMER / f"{model_name}_vs_n39_CANONICAL"
    tblout = prefix.with_suffix(".tblout")
    domtblout = prefix.with_suffix(".domtblout")
    stdout = prefix.with_suffix(".hmmsearch.txt")

    cmd = [
        "hmmsearch",
        "--cpu", "4",
        "--max",
        "--noali",
        "-E", REPORT_E,
        "--domE", REPORT_DOME,
        "--tblout", str(tblout),
        "--domtblout", str(domtblout),
        str(hmm_path),
        str(FASTA),
    ]

    print("\n== Running ==")
    print(model_name)
    print(" ".join(cmd))

    with stdout.open("w") as out:
        p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)

    if p.returncode != 0:
        fail(f"hmmsearch failed for {model_name}:\n{p.stderr}")

    seq_hits = parse_tblout(tblout)
    dom_hits = parse_domtblout(domtblout)

    pred_iE = {acc for acc, hit in dom_hits.items() if hit["independent_evalue"] <= PRIMARY_IE_THRESHOLD}
    pred_seqE = {acc for acc, hit in seq_hits.items() if hit["seq_evalue"] <= PRIMARY_IE_THRESHOLD}

    m_iE = compute_metrics(labels, pred_iE)
    m_seqE = compute_metrics(labels, pred_seqE)

    for criterion, metrics in [("best_domain_iE_le_1e-10", m_iE), ("sequence_E_le_1e-10", m_seqE)]:
        row = {
            "model": model_name,
            "criterion": criterion,
            "hmm_file": str(hmm_path),
            "hmm_sha256": sha256_file(hmm_path),
            "tblout": str(tblout),
            "domtblout": str(domtblout),
        }
        row.update(metrics)
        summary_rows.append(row)

    for acc, lab in labels.items():
        sh = seq_hits.get(acc)
        dh = dom_hits.get(acc)

        all_per_sequence.append({
            "model": model_name,
            "accession": acc,
            "truth": lab["binary_label"],
            "tfm_verdict": lab["tfm_verdict"],
            "length": lab["length"],
            "cys": lab["cys"],
            "domain_source": lab["domain_source"],
            "seq_evalue": sh["seq_evalue"] if sh else "",
            "seq_score": sh["seq_score"] if sh else "",
            "best_domain_iE": dh["independent_evalue"] if dh else "",
            "best_domain_cE": dh["conditional_evalue"] if dh else "",
            "domain_score": dh["domain_score"] if dh else "",
            "hmm_coords": f"{dh['hmm_from']}-{dh['hmm_to']}" if dh else "",
            "ali_coords": f"{dh['ali_from']}-{dh['ali_to']}" if dh else "",
            "pred_positive_iE_le_1e-10": acc in pred_iE,
            "pred_positive_seqE_le_1e-10": acc in pred_seqE,
        })

    commands[model_name] = {
        "command": " ".join(cmd),
        "tblout": str(tblout),
        "tblout_sha256": sha256_file(tblout),
        "domtblout": str(domtblout),
        "domtblout_sha256": sha256_file(domtblout),
        "stdout": str(stdout),
        "stdout_sha256": sha256_file(stdout),
    }

perseq_path = METRICS / "S8_HMMER_n39_CANONICAL_per_sequence.tsv"
summary_path = METRICS / "S8_HMMER_n39_CANONICAL_metrics.tsv"
manifest_path = METRICS / "S8_HMMER_n39_CANONICAL_manifest.json"

with perseq_path.open("w", newline="") as f:
    fieldnames = list(all_per_sequence[0].keys())
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(all_per_sequence)

with summary_path.open("w", newline="") as f:
    fieldnames = list(summary_rows[0].keys())
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    for r in summary_rows:
        w.writerow({k: fmt(v) for k, v in r.items()})

manifest = {
    "analysis": "S8 de novo HMM-only comparison on canonical n=39 mature-domain-compatible calibration set",
    "hmmer_version": hmmer_version,
    "primary_classification_rule": "HMM-positive if best independent-domain E-value <= 1e-10",
    "secondary_check": "sequence-level E-value <= 1e-10",
    "n39_fasta": str(FASTA),
    "n39_fasta_sha256": sha256_file(FASTA),
    "n39_labels": str(LABELS),
    "n39_labels_sha256": sha256_file(LABELS),
    "n_labels": len(labels),
    "label_counts": dict(Counter(v["binary_label"] for v in labels.values())),
    "models": {m: {"hmm": str(p), "sha256": sha256_file(p)} for m, p in HMM_MODELS.items()},
    "commands": commands,
    "outputs": {
        "per_sequence": str(perseq_path),
        "per_sequence_sha256": sha256_file(perseq_path),
        "metrics": str(summary_path),
        "metrics_sha256": sha256_file(summary_path),
    },
}
manifest_path.write_text(json.dumps(manifest, indent=2))

print("\n== Metrics ==")
print(summary_path.read_text())

print("== WROTE ==")
for p in [perseq_path, summary_path, manifest_path]:
    print(p, sha256_file(p))
