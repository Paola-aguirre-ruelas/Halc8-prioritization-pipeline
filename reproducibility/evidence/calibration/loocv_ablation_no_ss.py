"""
loocv_ablation_no_ss.py
════════════════════════════════════════════════════════════════════════════════
Secondary-structure ablation LOOCV for the HalC8 Detection Pipeline v2.2.

Scenario: Module 2 (secondary structure) completely removed.
  - p_ss   = 0.0 for every sequence (no SS points awarded)
  - ss_ok  = True for every sequence (HG2 treated as not assessable)
  - Module 1, Module 3, HG1, HG3, scoring weights: unchanged
  - Threshold optimisation: LOOCV-clean — same 32-combo grid search
    applied on 41 training sequences per fold, thresholds then used
    to classify the held-out sequence

Purpose:
  To quantify dependence of internal calibration performance on legacy
  GOR IV / Geneious 2021 secondary-structure labels. This is an end-to-end
  reproducibility sensitivity analysis, not a predictor-substitution test.

Input:   tabla5_tfm2021.csv  (same as loocv_validation.py)
Output:  ablation_no_ss_results.csv  +  printed comparison table

Usage:
    python3 loocv_ablation_no_ss.py

Run time: ~2 s on a modern laptop (32 combos × 42 folds = 1344 evaluations).
"""

import os, warnings
from datetime import datetime
import pandas as pd
import numpy as np

warnings.filterwarnings('ignore')

WORK_DIR = os.path.dirname(os.path.abspath(__file__))
TFM_CSV  = os.path.join(WORK_DIR, "tabla5_tfm2021.csv")
OUT_CSV  = os.path.join(WORK_DIR, "ablation_no_ss_results.csv")

MOTIF_PTS = {4: 3.0, 3: 2.0, 2: 1.0, 1: 0.5, 0: 0.0}

# Reference v2.2 reported metrics (for comparison table)
V22 = dict(tp=24, fp=1, tn=17, fn=0,
           sensitivity=1.000, specificity=0.944,
           precision=0.960, f1=0.980, kappa=0.951)

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")

# ── Helpers (identical to loocv_validation.py) ─────────────────────────────────

def normalize_verdict(v):
    v = str(v).strip().upper()
    if "MOST LIKELY" in v: return "MOST LIKELY"
    if "PROBABLE"   in v: return "PROBABLE"
    if "MARGINAL"   in v: return "MARGINAL"
    if "EXCLUDED"   in v: return "EXCLUDED"
    return "UNKNOWN"

def _parse_motif(motif_str):
    m = str(motif_str).strip().upper().replace('-','').replace('X','X')
    if not m or m in ('NAN','NOT FOUND',''): return 0
    s1 = 1 if m[0] in 'DENQ' else 0
    s2 = 1 if 'G' in m[:-1] else 0
    s3 = 1 if m.endswith('C') else 0
    s4 = 1 if (len(m) >= 2 and m[-2] in 'AVLI') else 0
    return s1 + s2 + s3 + s4

def score_row_no_ss(row):
    """
    Score one TFM row — Module 2 ablated.
    p_ss is forced to 0.0; ss_ok is forced to True (HG2 not assessable).
    Module 1 and Module 3 scoring unchanged.
    """
    def yn(col, w): return w if str(row.get(col,'NOT')).strip().upper()=='YES' else 0.0
    p_z1   = yn('hydrophob_first30', 1.5)
    p_z2   = yn('hydrophob_30_50',   1.5)
    p_z3   = yn('hydrophob_after50', 1.0)
    p_ss   = 0.0   # ABLATED: Module 2 zeroed
    ss_ok  = True  # ABLATED: HG2 not assessable
    ms     = _parse_motif(row.get('conserved_motif', ''))
    total  = round(p_z1 + p_z2 + p_z3 + p_ss + MOTIF_PTS[min(ms, 4)], 2)
    z1_yes = str(row.get('hydrophob_first30','NOT')).strip().upper() == 'YES'
    return total, ms, z1_yes, ss_ok

def predict(total, ms, z1_yes, ss_ok,
            z1_required=True, score_ml=6.5, score_prob=3.5):
    """
    Verdict rules — identical to v2.2.
    In ablation, ss_ok is always True so HG2 never triggers.
    """
    if z1_required and not z1_yes: return "EXCLUDED"   # HG1
    if not ss_ok:                  return "EXCLUDED"   # HG2 (never fires in ablation)
    if total >= score_ml  and ms == 4:              return "MOST LIKELY"
    if total >= score_prob and ms >= 2:             return "PROBABLE"
    if total >= 3.0 or (ms >= 3 and total >= 2.5): return "MARGINAL"
    return "EXCLUDED"

def cohen_kappa(actual, predicted):
    a = [1 if v in ('MOST LIKELY','PROBABLE') else 0 for v in actual]
    p = [1 if v in ('MOST LIKELY','PROBABLE') else 0 for v in predicted]
    n = len(a)
    po = sum(ai==pi for ai,pi in zip(a,p)) / n
    pe = (sum(a)/n)*(sum(p)/n) + ((n-sum(a))/n)*((n-sum(p))/n)
    return (po - pe) / (1 - pe) if pe < 1.0 else 0.0

# ── LOOCV with grid search ─────────────────────────────────────────────────────

def run_ablation_loocv(df):
    """
    LOOCV with Module 2 ablated.
    Threshold grid search is LOOCV-clean: optimised on 41 training sequences
    per fold, then applied to the held-out sequence.
    Grid is identical to loocv_validation.py:
      z1_required: [True, False]
      score_ml:    [6.5, 7.0, 7.5, 8.0]
      score_prob:  [3.5, 4.0, 4.5, 5.0]
    Note: with Module 2 ablated, max total score = 7.0 (M1=4.0 + M3=3.0),
    so score_ml = 7.5 and 8.0 are unreachable but included for symmetry.
    """
    Z1_OPTS    = [True, False]
    SCORE_ML   = [6.5, 7.0, 7.5, 8.0]
    SCORE_PROB = [3.5, 4.0, 4.5, 5.0]

    fold_results = []
    all_actuals, all_preds = [], []

    for i in range(len(df)):
        train = df[df.index != i]
        test  = df[df.index == i].iloc[0]
        actual = test['verdict_norm']

        # Grid search on 41 training sequences
        best_k, best_params = -99, (True, 6.5, 3.5)
        for z1r in Z1_OPTS:
            for sml in SCORE_ML:
                for spb in SCORE_PROB:
                    if spb >= sml: continue
                    preds = [predict(r.total, r.ms, r.z1_yes, r.ss_ok, z1r, sml, spb)
                             for r in train.itertuples()]
                    k = cohen_kappa(train['verdict_norm'].tolist(), preds)
                    if k > best_k:
                        best_k, best_params = k, (z1r, sml, spb)

        # Classify held-out sequence with fold-optimal thresholds
        z1r, sml, spb = best_params
        pred = predict(test.total, test.ms, test.z1_yes, test.ss_ok, z1r, sml, spb)

        all_actuals.append(actual)
        all_preds.append(pred)
        fold_results.append({
            'fold': i+1,
            'accession':        test.get('accession', f'seq_{i}'),
            'actual':           actual,
            'predicted':        pred,
            'match':            int(pred == actual),
            'best_z1_req':      z1r,
            'best_score_ml':    sml,
            'best_score_prob':  spb,
            'train_kappa':      round(best_k, 4),
            'total_score_no_ss': test.total,
            'motif_score':      test.ms,
        })

    return pd.DataFrame(fold_results), all_actuals, all_preds

# ── Metrics ────────────────────────────────────────────────────────────────────

def compute_metrics(all_actuals, all_preds):
    pos = {'MOST LIKELY', 'PROBABLE'}
    a = [1 if v in pos else 0 for v in all_actuals]
    p = [1 if v in pos else 0 for v in all_preds]
    n  = len(a)
    tp = sum(ai==1 and pi==1 for ai,pi in zip(a,p))
    tn = sum(ai==0 and pi==0 for ai,pi in zip(a,p))
    fp = sum(ai==0 and pi==1 for ai,pi in zip(a,p))
    fn = sum(ai==1 and pi==0 for ai,pi in zip(a,p))
    sens = tp/(tp+fn) if (tp+fn) else 0.0
    spec = tn/(tn+fp) if (tn+fp) else 0.0
    prec = tp/(tp+fp) if (tp+fp) else 0.0
    f1   = 2*prec*sens/(prec+sens) if (prec+sens) else 0.0
    kap  = cohen_kappa(all_actuals, all_preds)
    acc  = (tp+tn)/n
    kse  = np.sqrt(acc*(1-acc)/n)
    ci_lo = round(kap - 1.96*kse, 3)
    ci_hi = round(min(kap + 1.96*kse, 1.0), 3)
    return dict(n=n, tp=tp, tn=tn, fp=fp, fn=fn,
                sensitivity=sens, specificity=spec,
                precision=prec, f1=f1, kappa=kap,
                kappa_ci=(ci_lo, ci_hi))

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    log("=" * 65)
    log("HalC8 Pipeline v2.2 — No-SS Ablation LOOCV")
    log("Module 2 = 0.0 | HG2 = not assessable | thresholds LOOCV-clean")
    log("=" * 65)

    if not os.path.isfile(TFM_CSV):
        raise FileNotFoundError(f"Not found: {TFM_CSV}")

    df = pd.read_csv(TFM_CSV)
    df = df[df['tfm_verdict'].notna()].reset_index(drop=True)
    df = df[~df['accession'].str.contains('Annotated|reference', case=False, na=False)].reset_index(drop=True)
    log(f"Loaded {len(df)} sequences")

    df['verdict_norm'] = df['tfm_verdict'].apply(normalize_verdict)
    log(f"Verdict distribution: {dict(df['verdict_norm'].value_counts())}")

    # Pre-score all rows (ablated scoring, invariant across folds)
    df['total'], df['ms'], df['z1_yes'], df['ss_ok'] = zip(
        *df.apply(score_row_no_ss, axis=1))

    log("Running LOOCV (32 threshold combos × 42 folds)...")
    fold_df, all_actuals, all_preds = run_ablation_loocv(df)
    m = compute_metrics(all_actuals, all_preds)

    # ── Comparison table ───────────────────────────────────────────────────────
    def delta(key):
        diff = m[key] - V22[key]
        return f"{diff:+.3f}"

    log("")
    print(f"""
╔══════════════════════════════════════════════════════════════════╗
║       LOOCV comparison: v2.2 reported vs. No-SS ablation        ║
╠════════════════════╦══════════════════╦═════════════════════════╣
║ Metric             ║ v2.2 (reported)  ║ No-SS ablation (Δ)      ║
╠════════════════════╬══════════════════╬═════════════════════════╣
║ n                  ║ 42               ║ {m['n']:<24}║
║ TP                 ║ 24               ║ {m['tp']:<24}║
║ FP                 ║ 1                ║ {str(m['fp']) + '  ' + delta('fp'):<24}║
║ TN                 ║ 17               ║ {m['tn']:<24}║
║ FN                 ║ 0                ║ {str(m['fn']) + '  ' + delta('fn'):<24}║
║ Sensitivity        ║ 1.000            ║ {str(round(m['sensitivity'],3)) + '  ' + delta('sensitivity'):<24}║
║ Specificity        ║ 0.944            ║ {str(round(m['specificity'],3)) + '  ' + delta('specificity'):<24}║
║ Precision          ║ 0.960            ║ {str(round(m['precision'],3)) + '  ' + delta('precision'):<24}║
║ F1-score           ║ 0.980            ║ {str(round(m['f1'],3)) + '  ' + delta('f1'):<24}║
║ Cohen's κ          ║ 0.951            ║ {str(round(m['kappa'],3)) + '  ' + delta('kappa') + '  95% CI ' + str(m['kappa_ci']):<24}║
╚════════════════════╩══════════════════╩═════════════════════════╝
""")

    # ── Sequences whose verdict changed ────────────────────────────────────────
    changed = fold_df[fold_df['actual'] != fold_df['predicted']]
    log(f"Reclassified sequences: {len(changed)} of {len(fold_df)}")
    if len(changed):
        print("\nSequences reclassified vs. ground-truth label:")
        print(f"  {'Accession':<32} {'Actual':<14} {'Ablation pred':<14} {'Score(no-SS)'}")
        for _, row in changed.iterrows():
            print(f"  {row['accession']:<32} {row['actual']:<14} {row['predicted']:<14} {row['total_score_no_ss']:.1f}")
    else:
        print("  No reclassifications — ablation does not change any verdict.")

    # ── Interpretation hint ────────────────────────────────────────────────────
    delta_kappa = m['kappa'] - V22['kappa']
    delta_fp    = m['fp'] - V22['fp']
    print("\nInterpretation hint:")
    if abs(delta_kappa) < 0.05:
        print("  κ difference < 0.05 — ablation retained comparable discrimination.")
        print("  Suggested text: 'main classification behaviour driven primarily by")
        print("  hydrophobicity-zone architecture and motif geometry.'")
    elif delta_kappa < -0.05 and delta_fp > 0:
        print(f"  κ dropped {abs(delta_kappa):.3f}, FP increased by {delta_fp}.")
        print("  Suggested text: 'Module 2 contributes mainly to exclusion of")
        print("  N-terminally incompatible sequences; ablation reduces specificity")
        print("  while retaining sensitivity.'")
    else:
        print(f"  κ dropped {abs(delta_kappa):.3f}.")
        print("  Suggested text: 'legacy SS labels contributed materially to")
        print("  calibration-set discrimination; de novo candidates scored")
        print("  conservatively under NOT_PREDICTED policy.'")

    fold_df.to_csv(OUT_CSV, index=False)
    log(f"Per-sequence results saved: {OUT_CSV}")
    log("Done.")

if __name__ == "__main__":
    main()
