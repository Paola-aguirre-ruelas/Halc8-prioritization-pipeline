#!/usr/bin/env python3

from pathlib import Path
import math
import re
import pandas as pd


OUT = Path("audit_custom_hmm_tier2/results")
OUT.mkdir(parents=True, exist_ok=True)


def canon(x):
    s = str(x).strip().lstrip(">")
    m = re.search(r"(WP_\d+\.\d+)", s)
    if m:
        return m.group(1)
    return s.split()[0]


def parse_domtblout(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split(maxsplit=22)
            if len(parts) < 22:
                continue

            target = parts[0]
            tlen = int(parts[2])
            qlen = int(parts[5])
            full_eval = float(parts[6])
            full_score = float(parts[7])
            i_eval = float(parts[12])
            domain_bitscore = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])

            rows.append({
                "accession": canon(target),
                "raw_hmm_id": target,
                "target_len": tlen,
                "query_len": qlen,
                "full_evalue": full_eval,
                "full_score": full_score,
                "i_evalue": i_eval,
                "domain_bitscore": domain_bitscore,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "hmm_coverage": (hmm_to - hmm_from + 1) / qlen if qlen else math.nan,
                "alignment_length": ali_to - ali_from + 1,
            })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    return (
        df.sort_values(["accession", "i_evalue", "domain_bitscore"], ascending=[True, True, False])
          .groupby("accession", as_index=False)
          .first()
    )


def binary_from_verdict(v):
    s = str(v).strip().upper()
    s_norm = s.replace("-", "_").replace(" ", "_")

    positive = {
        "1",
        "TRUE",
        "YES",
        "POSITIVE",
        "SCORE_POSITIVE",
        "MOST_LIKELY",
        "PROBABLE",
    }

    negative = {
        "0",
        "FALSE",
        "NO",
        "NEGATIVE",
        "NOT_PRIORITIZED",
        "MARGINAL",
        "EXCLUDED",
    }

    if s_norm in positive:
        return 1
    if s_norm in negative:
        return 0

    if "MOST" in s and "LIKELY" in s:
        return 1
    if "PROBABLE" in s:
        return 1
    if "MARGINAL" in s:
        return 0
    if "EXCLUDED" in s:
        return 0
    if "NEGATIVE" in s:
        return 0
    if "POSITIVE" in s and "NOT" not in s:
        return 1

    return pd.NA


def metrics_at_threshold(y, x, thr):
    pred = (x >= thr).astype(int)
    y = y.astype(int)

    tp = int(((pred == 1) & (y == 1)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    n = tp + tn + fp + fn

    acc = (tp + tn) / n if n else math.nan
    sens = tp / (tp + fn) if (tp + fn) else math.nan
    spec = tn / (tn + fp) if (tn + fp) else math.nan
    prec = tp / (tp + fp) if (tp + fp) else math.nan
    f1 = 2 * prec * sens / (prec + sens) if prec == prec and sens == sens and (prec + sens) else math.nan

    pyes_true = (tp + fn) / n
    pyes_pred = (tp + fp) / n
    pno_true = (tn + fp) / n
    pno_pred = (tn + fn) / n
    pe = pyes_true * pyes_pred + pno_true * pno_pred
    kappa = (acc - pe) / (1 - pe) if (1 - pe) else math.nan

    denom = math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    mcc = ((tp*tn)-(fp*fn))/denom if denom else math.nan

    return {
        "threshold": thr,
        "n": n,
        "TP": tp,
        "TN": tn,
        "FP": fp,
        "FN": fn,
        "accuracy": acc,
        "sensitivity": sens,
        "specificity": spec,
        "precision": prec,
        "F1": f1,
        "kappa": kappa,
        "MCC": mcc,
    }


def auroc(y, x):
    d = pd.DataFrame({"y": y, "x": x}).dropna()
    if d["y"].nunique() < 2:
        return math.nan
    pos = d[d.y == 1]
    neg = d[d.y == 0]
    ranks = d["x"].rank(method="average")
    sum_pos = ranks[d.y == 1].sum()
    return float((sum_pos - len(pos)*(len(pos)+1)/2) / (len(pos)*len(neg)))


def threshold_scan(df, ycol, scorecol):
    d = df[[ycol, scorecol]].dropna()
    rows = [metrics_at_threshold(d[ycol], d[scorecol], thr) for thr in sorted(d[scorecol].unique())]
    out = pd.DataFrame(rows)
    return out.sort_values(["kappa", "F1", "MCC"], ascending=False)


def analyze_calibration():
    hmm = parse_domtblout("audit_custom_hmm_tier2/hmmer/calibration_42_vs_HalC8_custom_v2.2.domtblout")
    labels = pd.read_csv("./calibration_42_labels.csv")
    tier2 = pd.read_csv("./results/modular_ablation_per_sequence.csv")

    labels["accession"] = labels["accession"].map(canon)
    tier2["accession"] = tier2["accession"].map(canon)

    labels["reference_binary"] = labels["binary_label"].map(binary_from_verdict)
    
    if labels["reference_binary"].isna().any():
    	bad = labels.loc[labels["reference_binary"].isna(), ["accession", "binary_label"]]
    	raise SystemExit(
            "Unmapped calibration binary_label values:\n" + bad.to_string(index=False)
    	)

    labels["reference_binary"] = labels["reference_binary"].astype(int)

    # Full Tier 2 score and no-SS score from reconstructed audit.
    keep = [
        "accession",
        "actual_binary",
        "full_score_calc",
        "no_ss_score_calc",
        "pred_full_reconstructed",
        "pred_no_ss_reconstructed",
        "M1_hydro",
        "M2_structure",
        "M3_motif",
        "hg1_pass",
        "hg2_pass",
        "motif_score_4_calc",
        "genus",
        "tier",
    ]
    keep = [c for c in keep if c in tier2.columns]

    df = labels.merge(tier2[keep], on="accession", how="left")
    df = df.merge(hmm, on="accession", how="left")

    df.to_csv(OUT / "calibration_42_custom_hmm_vs_tier2.tsv", sep="\t", index=False)

    scan = threshold_scan(df, "reference_binary", "domain_bitscore")
    scan.to_csv(OUT / "calibration_42_hmm_bitscore_threshold_scan.tsv", sep="\t", index=False)

    best = scan.iloc[0].to_dict()
    min_pos = df[df.reference_binary == 1]["domain_bitscore"].min()
    max_neg = df[df.reference_binary == 0]["domain_bitscore"].max()

    spearman_full = df["domain_bitscore"].rank().corr(df["full_score_calc"].rank())
    spearman_noss = df["domain_bitscore"].rank().corr(df["no_ss_score_calc"].rank())

    summary = {
        "dataset": "calibration_42",
        "n": len(df),
        "hmm_rows_merged": int(df["domain_bitscore"].notna().sum()),
        "AUROC_HMM_bitscore_vs_reference_binary": auroc(df["reference_binary"], df["domain_bitscore"]),
        "Spearman_HMM_bitscore_vs_full_Tier2_score": spearman_full,
        "Spearman_HMM_bitscore_vs_noSS_Tier2_score": spearman_noss,
        "min_positive_domain_bitscore": min_pos,
        "max_negative_domain_bitscore": max_neg,
        "HMM_bitscore_perfectly_separates_reference_binary": bool(min_pos > max_neg),
        **{f"best_HMM_threshold_{k}": v for k, v in best.items()},
    }

    pd.DataFrame([summary]).to_csv(OUT / "calibration_42_summary.tsv", sep="\t", index=False)


def analyze_refseq():
    hmm = parse_domtblout("audit_custom_hmm_tier2/hmmer/refseq_18_vs_HalC8_custom_v2.2.domtblout")
    tier2 = pd.read_csv("./halobacteria_18_pipeline_scores.csv")

    tier2["accession_base"] = tier2["accession"].map(canon)
    tier2["scorer_binary"] = tier2["verdict_v2"].map(binary_from_verdict).astype(int)

    df = tier2.merge(hmm, left_on="accession_base", right_on="accession", how="left", suffixes=("_tier2", "_hmm"))

    df.to_csv(OUT / "refseq_18_custom_hmm_vs_tier2.tsv", sep="\t", index=False)

    scan = threshold_scan(df, "scorer_binary", "domain_bitscore")
    scan.to_csv(OUT / "refseq_18_hmm_bitscore_threshold_scan.tsv", sep="\t", index=False)

    best = scan.iloc[0].to_dict()
    min_pos = df[df.scorer_binary == 1]["domain_bitscore"].min()
    max_neg = df[df.scorer_binary == 0]["domain_bitscore"].max()

    spearman = df["domain_bitscore"].rank().corr(df["total_score_9"].rank())

    ranked = df.sort_values("domain_bitscore", ascending=False)
    # Figure 3C annotation: Halorussus identity in the frozen RefSeq deployment set.
    # This affects marker shape only; it does not alter HMM, Tier 2 score, gates, or metrics.
    halorussus_accessions = {
        "WP_276279820.1",
        "WP_253478286.1",
        "WP_166035560.1",
    }
    if "accession_base" not in ranked.columns:
        raise KeyError("Expected column 'accession_base' not found in ranked RefSeq table.")
    ranked["is_halorussus"] = ranked["accession_base"].apply(
        lambda x: "TRUE" if str(x).strip() in halorussus_accessions else "FALSE"
    )

    ranked.to_csv(OUT / "refseq_18_ranked_by_custom_hmm_bitscore.tsv", sep="\t", index=False)

    summary = {
        "dataset": "refseq_18",
        "n": len(df),
        "hmm_rows_merged": int(df["domain_bitscore"].notna().sum()),
        "n_score_positive": int(df["scorer_binary"].sum()),
        "n_not_prioritized": int((df["scorer_binary"] == 0).sum()),
        "AUROC_HMM_bitscore_vs_Tier2_binary": auroc(df["scorer_binary"], df["domain_bitscore"]),
        "Spearman_HMM_bitscore_vs_Tier2_score": spearman,
        "min_score_positive_domain_bitscore": min_pos,
        "max_not_prioritized_domain_bitscore": max_neg,
        "HMM_bitscore_perfectly_separates_Tier2_binary": bool(min_pos > max_neg),
        **{f"best_HMM_threshold_{k}": v for k, v in best.items()},
    }

    pd.DataFrame([summary]).to_csv(OUT / "refseq_18_summary.tsv", sep="\t", index=False)


def main():
    analyze_calibration()
    analyze_refseq()

    print("Wrote:")
    for f in sorted(OUT.glob("*.tsv")):
        print(f)


if __name__ == "__main__":
    main()
