"""Check core artifacts and reproduce three archived scoring tables."""
from pathlib import Path
import hashlib
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "reproducibility"
for line in (R / "SHA256SUMS.txt").read_text().splitlines():
    digest, name = line.split("  ", 1)
    if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
        raise SystemExit(f"Hash mismatch: {name}")
print("Artifact integrity: PASS", flush=True)

CASES = [
    ("refseq18_for_scoring.csv", "none", "halobacteria_18_pipeline_scores_CLEAN.csv"),
    ("calibration42_mature_domains.csv", "none", "calibration42_mature_NONE_scores.csv"),
    ("calibration42_mature_domains.csv", "chou_fasman", "calibration42_mature_CF_scores.csv"),
]
with tempfile.TemporaryDirectory() as directory:
    for input_name, mode, expected in CASES:
        output = Path(directory) / expected
        result = subprocess.run([sys.executable, str(ROOT / "halc8_score_domains.py"),
            "--input", str(R / "data" / input_name), "--output", str(output),
            "--input-type", "domain", "--ss-mode", mode], capture_output=True, text=True)
        if result.returncode:
            raise SystemExit(result.stdout + result.stderr)
        if output.read_bytes() != (R / "expected" / expected).read_bytes():
            raise SystemExit(f"CSV mismatch: {expected}")
        print(f"PASS: {expected}")
