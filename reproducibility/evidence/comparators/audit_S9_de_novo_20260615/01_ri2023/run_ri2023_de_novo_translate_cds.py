from pathlib import Path
import csv, json, hashlib, subprocess, sys

ROOT = Path("audit_S9_de_novo_20260615")
OUT = ROOT / "01_ri2023"
OUT.mkdir(parents=True, exist_ok=True)

NT_FASTA = Path("OP548098.1.fasta")
HMM = Path("HalC8_custom_v2.2.hmm")
SCORER = Path("halc8_score_domains.py")

GENETIC_CODE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

EXPECTED_PROTEIN_FROM_OLD_TABLE = (
    "DIDITGCSACKYAAGQVCTIGCSAAGGFICGLLGITIPVAGLSCLGFVEIVCTVADEYSGCGDAVAKEACNRVGLC"
)

def fail(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

def sha256_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def sha256_text(s):
    return hashlib.sha256(s.encode()).hexdigest()

def read_single_fasta(path):
    header = None
    seq = []
    for line in path.read_text(errors="ignore").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                fail("FASTA contains more than one record")
            header = line[1:]
        else:
            seq.append(line)
    if header is None:
        fail("FASTA has no header")
    sequence = "".join(seq).replace(" ", "").replace("\t", "").upper().replace("U", "T")
    return header, sequence

def translate_cds(nt):
    if set(nt) - set("ACGT"):
        fail(f"Non-DNA characters found: {sorted(set(nt)-set('ACGT'))}")
    if len(nt) % 3 != 0:
        fail(f"CDS length is not divisible by 3: {len(nt)}")
    aas = []
    for i in range(0, len(nt), 3):
        codon = nt[i:i+3]
        aa = GENETIC_CODE.get(codon)
        if aa is None:
            fail(f"Unknown codon: {codon}")
        aas.append(aa)
    prot = "".join(aas)
    if "*" in prot[:-1]:
        fail(f"Internal stop codon in translated protein: {prot}")
    if prot.endswith("*"):
        prot = prot[:-1]
    return prot

def parse_best_domtblout(path):
    hits = []
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            hits.append({
                "target": parts[0],
                "query": parts[3],
                "full_evalue": float(parts[6]),
                "full_score": float(parts[7]),
                "conditional_evalue": float(parts[11]),
                "independent_evalue": float(parts[12]),
                "domain_score": float(parts[13]),
                "hmm_from": int(parts[15]),
                "hmm_to": int(parts[16]),
                "ali_from": int(parts[17]),
                "ali_to": int(parts[18]),
                "env_from": int(parts[19]),
                "env_to": int(parts[20]),
                "acc": float(parts[21]),
            })
    if not hits:
        fail("No HMMER domain hit found after CDS translation")
    return sorted(hits, key=lambda x: (x["independent_evalue"], -x["domain_score"]))[0]

for p in [NT_FASTA, HMM, SCORER]:
    if not p.exists():
        fail(f"Missing required file: {p}")

header, nt_seq = read_single_fasta(NT_FASTA)
aa_seq = translate_cds(nt_seq)
acc = "OP548098.1"

protein_fasta = OUT / "ri2023_OP548098_translated_CDS_protein.faa"
input_csv = OUT / "ri2023_OP548098_translated_CDS_CANONICAL_input.csv"
scores_csv = OUT / "ri2023_OP548098_translated_CDS_CANONICAL_scores.csv"
tblout = OUT / "ri2023_OP548098_translated_CDS_vs_HalC8_custom.tblout"
domtblout = OUT / "ri2023_OP548098_translated_CDS_vs_HalC8_custom.domtblout"
hmmer_txt = OUT / "ri2023_OP548098_translated_CDS_vs_HalC8_custom.hmmsearch.txt"
summary_tsv = OUT / "ri2023_OP548098_translated_CDS_CANONICAL_summary.tsv"
manifest_json = OUT / "ri2023_OP548098_translated_CDS_CANONICAL_manifest.json"

with protein_fasta.open("w") as f:
    f.write(
        f">{acc}|translated_CDS|len{len(aa_seq)}|cys{aa_seq.count('C')}|"
        f"ntsha256:{sha256_text(nt_seq)}|aasha256:{sha256_text(aa_seq)}\n"
    )
    f.write(aa_seq + "\n")

with input_csv.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["accession", "sequence", "category", "notes"])
    w.writeheader()
    w.writerow({
        "accession": acc,
        "sequence": aa_seq,
        "category": "post_calibration_positive_Ri2023_translated_CDS",
        "notes": header,
    })

hmmer_cmd = [
    "hmmsearch",
    "--cpu", "4",
    "--max",
    "--noali",
    "-E", "1000",
    "--domE", "1000",
    "--tblout", str(tblout),
    "--domtblout", str(domtblout),
    str(HMM),
    str(protein_fasta),
]

with hmmer_txt.open("w") as out:
    p = subprocess.run(hmmer_cmd, stdout=out, stderr=subprocess.PIPE, text=True)
if p.returncode != 0:
    fail(f"HMMER failed:\n{p.stderr}")

score_cmd = [
    "python3", str(SCORER),
    "--input", str(input_csv),
    "--output", str(scores_csv),
    "--input-type", "domain",
    "--ss-mode", "none",
]

p = subprocess.run(score_cmd, capture_output=True, text=True)
(OUT / "ri2023_OP548098_translated_CDS_CANONICAL_score_stdout.txt").write_text(p.stdout)
(OUT / "ri2023_OP548098_translated_CDS_CANONICAL_score_stderr.txt").write_text(p.stderr)
if p.returncode != 0:
    fail(f"Scorer failed:\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")

best = parse_best_domtblout(domtblout)

with scores_csv.open() as f:
    rows = list(csv.DictReader(f))
if len(rows) != 1:
    fail(f"Expected one scorer row, observed {len(rows)}")
score = rows[0]

summary = {
    "source": "Ri et al. 2023 post-calibration positive audit",
    "accession": acc,
    "original_fasta_header": header,
    "input_unit": "translated nucleotide CDS; translated in-frame from OP548098.1.fasta",
    "nt_length": len(nt_seq),
    "aa_length": len(aa_seq),
    "aa_cys": aa_seq.count("C"),
    "nt_sha256": sha256_text(nt_seq),
    "aa_sha256": sha256_text(aa_seq),
    "translated_protein": aa_seq,
    "matches_old_table_sequence": aa_seq == EXPECTED_PROTEIN_FROM_OLD_TABLE,
    "hmm": str(HMM),
    "hmm_sha256": sha256_file(HMM),
    "hmmer_independent_evalue": best["independent_evalue"],
    "hmmer_domain_score": best["domain_score"],
    "hmmer_hmm_coords": f"{best['hmm_from']}-{best['hmm_to']}",
    "hmmer_alignment_coords": f"{best['ali_from']}-{best['ali_to']}",
    "scorer": str(SCORER),
    "scorer_sha256": sha256_file(SCORER),
    "input_type": score["input_type"],
    "ss_mode": "none",
    "z1_mean": score["z1_mean"],
    "z2_mean": score["z2_mean"],
    "z3_mean": score["z3_mean"],
    "pts_hydro_z1": score["pts_hydro_z1"],
    "pts_hydro_z2": score["pts_hydro_z2"],
    "pts_hydro_z3": score["pts_hydro_z3"],
    "pts_structure": score["pts_structure"],
    "motif_score_4": score["motif_score_4"],
    "pts_motif": score["pts_motif"],
    "total_score_9": score["total_score_9"],
    "verdict_v2": score["verdict_v2"],
    "hg_z1_pass": score["hg_z1_pass"],
    "hg_ss_nterm_ok": score["hg_ss_nterm_ok"],
    "hg_motif_anchor": score["hg_motif_anchor"],
    "conserved_motif": score["conserved_motif"],
}

with summary_tsv.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary.keys()), delimiter="\t")
    w.writeheader()
    w.writerow(summary)

manifest = {
    "analysis": "S9 Part A de novo Ri2023 post-calibration positive audit using translated OP548098.1 CDS",
    "input_nt_fasta": str(NT_FASTA),
    "input_nt_fasta_sha256": sha256_file(NT_FASTA),
    "nt_sequence_sha256": sha256_text(nt_seq),
    "translated_protein_fasta": str(protein_fasta),
    "translated_protein_fasta_sha256": sha256_file(protein_fasta),
    "aa_sequence_sha256": sha256_text(aa_seq),
    "hmm": str(HMM),
    "hmm_sha256": sha256_file(HMM),
    "scorer": str(SCORER),
    "scorer_sha256": sha256_file(SCORER),
    "hmmer_command": " ".join(hmmer_cmd),
    "score_command": " ".join(score_cmd),
    "outputs": {
        "input_csv": str(input_csv),
        "scores_csv": str(scores_csv),
        "tblout": str(tblout),
        "domtblout": str(domtblout),
        "hmmer_txt": str(hmmer_txt),
        "summary_tsv": str(summary_tsv),
    },
}

manifest_json.write_text(json.dumps(manifest, indent=2))

print("== Ri2023 de novo translated-CDS summary ==")
for k in [
    "accession",
    "nt_length",
    "aa_length",
    "aa_cys",
    "matches_old_table_sequence",
    "hmmer_independent_evalue",
    "hmmer_domain_score",
    "hmmer_hmm_coords",
    "hmmer_alignment_coords",
    "total_score_9",
    "motif_score_4",
    "conserved_motif",
    "verdict_v2",
]:
    print(f"{k}: {summary[k]}")

print("\n== Translated protein ==")
print(aa_seq)

print("\n== WROTE ==")
for p in [protein_fasta, input_csv, scores_csv, tblout, domtblout, hmmer_txt, summary_tsv, manifest_json]:
    print(p, sha256_file(p))
