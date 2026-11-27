#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_figure2.py
═══════════════════════════════════════════════════════════════════════════════
Figure 2 — Internal calibration, threshold stability, no-SS ablation,
and custom-HMM-only comparator for HalC8-Prioritizer.

Designed from scratch for BMC Bioinformatics submission.
Same navy + teal + neutral grey palette as Figure 1, journal-grade and
colour-blind safe. All plotted numeric results are read from archived TSVs. No synthetic κ grid is used.

Layout: A calibration-frame agreement; B operating-point stability and calibration-frame no-SS ablation; C HMM bit score versus Tier 2 score.
INPUTS (default: ./inputs/)
  figure2_panelA_metrics.tsv
  figure2_panelB_summary.tsv
  figure2_panelB_no_ss_ablation.tsv (Supplementary Table S9 Panel C)
  figure2_panelC_hmm_vs_tier2.tsv
  calibration_42_summary.tsv and refseq_18_summary.tsv (archived HMM audit)

OUTPUTS  (default: current Figure2_calibration_hmm_audit_source/ directory)
  Figure_2.svg   — fonts as live text (editable in Illustrator/Inkscape)
  Figure_2.pdf   — vectorial, BMC-grade

USAGE
  python3 generate_figure2.py
  python3 generate_figure2.py --indir inputs --outdir .
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle


# ──────────────────────────────────────────────────────────────────────────────
# Palette  —  uniform with Figure 1 (BMC Bioinformatics compatible)
# ──────────────────────────────────────────────────────────────────────────────
NAVY_900   = "#0d2240"
NAVY_700   = "#1a3a5c"
NAVY_500   = "#3d5a80"
NAVY_300   = "#8aa9d0"
NAVY_100   = "#dfe8f4"
NAVY_050   = "#f0f5fb"

TEAL_700   = "#1f6f6f"
TEAL_500   = "#2a8b8b"
TEAL_300   = "#7bbcbc"
TEAL_100   = "#dff0f0"

GREY_700   = "#3a3a3a"
GREY_500   = "#6b7280"
GREY_300   = "#aab1ba"
GREY_100   = "#f4f5f7"

INK        = "#1a1a1a"
WHITE      = "#ffffff"


# ──────────────────────────────────────────────────────────────────────────────
# Matplotlib styling
# ──────────────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":        9.5,
    "axes.titlesize":   11,
    "axes.labelsize":   9.5,
    "xtick.labelsize":  8.5,
    "ytick.labelsize":  8.5,
    "legend.fontsize":  8.5,
    "figure.dpi":       150,
    # Editable text in PDF and SVG
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "svg.fonttype":     "none",
    "axes.spines.top":  False,
    "axes.spines.right": False,
})

# Sequential κ heatmap colormap — pale → teal → navy
KAPPA_CMAP = LinearSegmentedColormap.from_list(
    "navy_teal",
    [(0.0,  NAVY_050),
     (0.25, TEAL_100),
     (0.55, TEAL_500),
     (0.85, TEAL_700),
     (1.0,  NAVY_700)],
    N=256,
)


# ──────────────────────────────────────────────────────────────────────────────
# I/O
# ──────────────────────────────────────────────────────────────────────────────
def _read_wide_tsv(path: Path) -> dict:
    """Read a wide TSV with a single data row → lower-cased key dict."""
    if not path.is_file():
        sys.exit(f"ERROR: required input not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if df.empty:
        sys.exit(f"ERROR: {path.name} has no data rows.")
    if len(df) > 1:
        print(f"  [warn] {path.name} has {len(df)} rows; using the first.",
              file=sys.stderr)
    row = df.iloc[0]
    return {str(c).strip().lower(): str(row[c]).strip() for c in df.columns}


def _strict_int(v):
    try: return int(round(float(v)))
    except Exception: return None


def _strict_float(v):
    try: return float(v)
    except Exception: return None


def _strict_bool(v):
    s = str(v).strip().lower()
    if s in {"true", "t", "yes", "y", "1"}:  return True
    if s in {"false", "f", "no", "n", "0"}:  return False
    return None


def _require(d: dict, keys: list, where: str):
    missing = [k for k in keys if k not in d]
    if missing:
        sys.exit(f"ERROR: {where} missing required keys: {missing}\n"
                 f"       Found columns: {sorted(d.keys())}")


def load_panel_a(path: Path) -> dict:
    d = _read_wide_tsv(path)
    needed = ["n", "tp", "tn", "fp", "fn",
              "sensitivity", "specificity", "precision", "f1", "kappa"]
    _require(d, needed, path.name)
    return {
        "n":           _strict_int(d["n"]),
        "tp":          _strict_int(d["tp"]),
        "tn":          _strict_int(d["tn"]),
        "fp":          _strict_int(d["fp"]),
        "fn":          _strict_int(d["fn"]),
        "kappa":       _strict_float(d["kappa"]),
        "sensitivity": _strict_float(d["sensitivity"]),
        "specificity": _strict_float(d["specificity"]),
        "precision":   _strict_float(d["precision"]),
        "f1":          _strict_float(d["f1"]),
    }


def load_panel_b_summary(path: Path) -> dict:
    d = _read_wide_tsv(path)
    # Accept "best_*" aliases produced by the canonical LOOCV pipeline.
    aliases = {
        "deployed_score_ml":   ("best_score_ml",),
        "deployed_score_prob": ("best_score_prob",),
        "z1_required":         ("best_z1_req",),
        "folds_evaluated":     (),
        "folds_same_operating_point": (),
        "deployed_kappa":      (),
    }
    canon = {}
    for k, alts in aliases.items():
        if k in d:
            canon[k] = d[k]
            continue
        for alt in alts:
            if alt in d:
                canon[k] = d[alt]
                break
    _require(canon, list(aliases.keys()), path.name)
    z1 = _strict_bool(canon["z1_required"])
    if z1 is None:
        z1 = True
    return {
        "ml":        _strict_float(canon["deployed_score_ml"]),
        "pr":        _strict_float(canon["deployed_score_prob"]),
        "z1":        z1,
        "n_folds":   _strict_int(canon["folds_evaluated"]),
        "n_agree":   _strict_int(canon["folds_same_operating_point"]),
        "kappa":     _strict_float(canon["deployed_kappa"]),
        # train κ range — used to detect anomalies, not rendered
        "kappa_min": _strict_float(d.get("train_kappa_min", "")) if "train_kappa_min" in d else None,
        "kappa_max": _strict_float(d.get("train_kappa_max", "")) if "train_kappa_max" in d else None,
    }


def load_panel_b_ablation(path: Path) -> dict:
    df = pd.read_csv(path, sep="\t", comment="#")
    required = {"metric", "default_mode", "no_ss_mode", "delta"}
    if not required.issubset(df.columns) or df["metric"].duplicated().any():
        raise ValueError(f"Invalid ablation input: {path}")
    df = df.set_index("metric")
    expected = ["Cohen's kappa", "Sensitivity", "Specificity", "Precision", "F1"]
    if set(df.index) != set(expected):
        raise ValueError(f"Unexpected ablation metrics: {path}")
    return {k: tuple(float(df.loc[k, col]) for col in ("default_mode", "no_ss_mode", "delta")) for k in expected}


def load_panel_c(path: Path) -> pd.DataFrame:
    """
    Load Panel C scatter data: one row per RefSeq sequence-novel domain.

    Required columns:
      accession, hmm_bitscore, tier2_score, verdict, is_score_positive
    Optional:
      is_halorussus, label
    """
    if not path.is_file():
        sys.exit(f"ERROR: required input not found: {path}")
    df = pd.read_csv(path, sep="\t", comment="#")
    df.columns = [c.strip() for c in df.columns]
    required = {"accession", "hmm_bitscore", "tier2_score",
                "verdict", "is_score_positive"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: {path.name} missing required columns: {sorted(missing)}")
    df["hmm_bitscore"]      = df["hmm_bitscore"].astype(float)
    df["tier2_score"]       = df["tier2_score"].astype(float)
    df["is_score_positive"] = df["is_score_positive"].astype(int)
    if "is_halorussus" not in df.columns:
        df["is_halorussus"] = False
    else:
        df["is_halorussus"] = (df["is_halorussus"].astype(str).str.strip().str.lower()
                               .isin(["true", "t", "yes", "y", "1"]))
    if "label" not in df.columns:
        df["label"] = ""
    df["label"] = df["label"].fillna("").astype(str)
    return df.reset_index(drop=True)


# ──────────────────────────────────────────────────────────────────────────────
# Panel A — Internal binary agreement
# ──────────────────────────────────────────────────────────────────────────────
def panel_a(ax_cm, ax_box, a: dict):
    """Confusion matrix (top) + metrics card (bottom)."""

    # Cell layout (canonical orientation, legacy label rows × pipeline call cols)
    # [TN  FP]
    # [FN  TP]
    cm = np.array([[a["tn"], a["fp"]],
                   [a["fn"], a["tp"]]], dtype=float)
    fill = np.array([[1.0, 0.15],
                     [0.15, 1.0]])
    cmap = LinearSegmentedColormap.from_list(
        "navy_two", [NAVY_100, NAVY_700], N=256)
    ax_cm.imshow(fill, cmap=cmap, vmin=0, vmax=1.05, aspect="equal")

    for i in range(2):
        for j in range(2):
            color = WHITE if fill[i, j] > 0.5 else NAVY_900
            ax_cm.text(j, i, str(int(cm[i, j])),
                       ha="center", va="center",
                       fontsize=26, fontweight="bold", color=color)

    # Short labels in the figure; full definition lives in the caption.
    neg_lab = "Not prioritized"
    pos_lab = "Score-positive"
    ax_cm.set_xticks([0, 1])
    ax_cm.set_yticks([0, 1])
    ax_cm.set_xticklabels([neg_lab, pos_lab], fontsize=9)
    ax_cm.set_yticklabels([neg_lab, pos_lab], fontsize=9,
                          rotation=90, va="center")
    ax_cm.tick_params(axis="x", which="both", length=0, pad=6)
    ax_cm.tick_params(axis="y", which="both", length=0, pad=20)
    ax_cm.set_xlabel("Pipeline binary call",
                     fontsize=9.5, fontweight="bold", labelpad=8,
                     color=NAVY_700)
    ax_cm.set_ylabel("Prior in-silico reference label",
                     fontsize=9.5, fontweight="bold", labelpad=10,
                     color=NAVY_700)
    # Xticks at top (class labels above the matrix), but xlabel at BOTTOM
    # so it never collides with the title up top.
    ax_cm.xaxis.tick_top()
    for s in ("top", "right", "bottom", "left"):
        ax_cm.spines[s].set_visible(False)

    # Title shifted left in the axes via explicit x; y matches Panel B.
    ax_cm.set_title("A   Calibration-frame agreement",
                    fontsize=12, fontweight="bold", loc="left",
                    pad=0, x=-0.13, y=1.20, color=NAVY_700)
    # Anchor the matrix at the SOUTH of its GridSpec cell so it sits closer
    # to the metrics card below (matches the user's "matrix lower" request).
    ax_cm.set_anchor("S")

    # Metrics card (below)
    ax_box.set_axis_off()
    lines = [
        ("Calibration",                 f"n = {a['n']}"),
        ("TP / TN / FP / FN",           f"{a['tp']} / {a['tn']} / {a['fp']} / {a['fn']}"),
        ("Cohen's κ",                   f"{a['kappa']:.3f}"),
        ("Sensitivity",                 f"{a['sensitivity']:.3f}"),
        ("Specificity",                 f"{a['specificity']:.3f}"),
        ("Precision",                   f"{a['precision']:.3f}"),
        ("F1",                          f"{a['f1']:.3f}"),
    ]
    # Card pinned high in its axes so it sits close to the matrix above.
    _two_col_card(ax_box, lines,
                  x=0.02, y_top=1.10, width=0.96, line_h=0.115,
                  font_size=8.8)
    # Italic caveat removed per spec; this information lives in the caption.


def _two_col_card(ax, pairs, x=0.02, y_top=0.96, width=0.96,
                  line_h=0.115, font_size=8.6,
                  facecolor=NAVY_050, edgecolor=NAVY_300, label_frac=0.55):
    """Two-column card with dotted leaders. Shared layout component."""
    from matplotlib.lines import Line2D
    n = len(pairs)
    card_h = line_h * n + 0.06
    card_y = y_top - card_h
    ax.add_patch(FancyBboxPatch(
        (x, card_y), width, card_h,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        facecolor=facecolor, edgecolor=edgecolor, linewidth=0.9,
        transform=ax.transAxes, clip_on=False, zorder=1))
    pad_l, pad_r = 0.022, 0.022
    x_label  = x + pad_l
    x_value  = x + width - pad_r
    leader_left  = x_label + width * (label_frac - pad_l - 0.02)
    leader_right = x_value - 0.04
    for i, (label, value) in enumerate(pairs):
        y = y_top - 0.04 - (i + 0.5) * line_h
        ax.text(x_label, y, str(label),
                transform=ax.transAxes, ha="left", va="center",
                fontsize=font_size, color=INK, zorder=3)
        ax.text(x_value, y, str(value),
                transform=ax.transAxes, ha="right", va="center",
                fontsize=font_size, color=INK,
                family="monospace", fontweight="bold", zorder=3)
        leader = Line2D([leader_left, leader_right], [y, y],
                        transform=ax.transAxes,
                        color=GREY_300, linewidth=0.7,
                        linestyle=(0, (1.2, 2.2)),
                        alpha=0.7, zorder=2, clip_on=False)
        ax.add_line(leader)


# ──────────────────────────────────────────────────────────────────────────────
# Panel B — Threshold stability and calibration-frame no-SS ablation
# ──────────────────────────────────────────────────────────────────────────────
def panel_b(ax, ax_box, b_summary: dict, ablation: dict):
    """Threshold stability + no-SS ablation cost.

    Replaces the previous κ heatmap with a paired-bar comparison of the
    default scoring mode vs the no-SS ablation. Three metrics shown:
    Cohen's κ, Specificity, Sensitivity (sensitivity included to show that
    the ablation cost is specificity, not recall).
    """
    metrics = [(label, *ablation[key][:2]) for label, key in
               [("Cohen's κ", "Cohen's kappa"), ("Specificity", "Specificity"), ("Sensitivity", "Sensitivity")]]
    labels      = [m[0] for m in metrics]
    full_vals   = [m[1] for m in metrics]
    noss_vals   = [m[2] for m in metrics]

    y_pos    = np.arange(len(metrics))
    bar_h    = 0.36

    # Horizontal paired bars (full = navy, no-SS = teal)
    ax.barh(y_pos + bar_h/2, full_vals,
            height=bar_h, color=NAVY_700, edgecolor=WHITE, linewidth=0.6,
            label="Full calibration", zorder=3)
    ax.barh(y_pos - bar_h/2, noss_vals,
            height=bar_h, color=TEAL_500, edgecolor=WHITE, linewidth=0.6,
            label="No-SS ablation", zorder=3)

    # Annotate bar values
    for i, (lab, full_v, noss_v) in enumerate(metrics):
        ax.text(full_v + 0.012, i + bar_h/2, f"{full_v:.3f}",
                va="center", ha="left", fontsize=8.4,
                fontweight="bold", color=NAVY_700)
        ax.text(noss_v + 0.012, i - bar_h/2, f"{noss_v:.3f}",
                va="center", ha="left", fontsize=8.4,
                fontweight="bold", color=TEAL_700)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9.3)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.18)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Metric value", fontsize=9.3, fontweight="bold",
                  color=NAVY_700)
    ax.grid(True, axis="x", alpha=0.25, linewidth=0.5, linestyle=":")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # Title y matches Panel A's so both align at the same absolute height.
    ax.set_title("B   Threshold stability",
                 fontsize=12, fontweight="bold", loc="left",
                 pad=0, y=1.20, color=NAVY_700)
    # Subtitle removed per spec; n = 42 calibration set context goes in caption.

    # Legend ABOVE the chart, lowered slightly to leave more breathing room
    # between it and the title.
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.005),
              ncol=2, fontsize=8, frameon=True,
              framealpha=0.95, edgecolor=NAVY_300, borderpad=0.4,
              columnspacing=1.5, handletextpad=0.5)

    # ── Summary card (ax_box) ─────────────────────────────────────────
    ax_box.set_axis_off()
    op_ml = b_summary["ml"]
    op_pr = b_summary["pr"]
    pairs = [
        ("LOOCV folds agreeing",       f"{b_summary['n_agree']} / {b_summary['n_folds']}"),
        ("Fixed score_ML",             f"{op_ml:.1f}"),
        ("Fixed score_PROB",           f"{op_pr:.1f}"),
        ("HG1 gate",                   "enabled" if b_summary['z1'] else "disabled"),
        ("Δκ on SS ablation",          ("{:+.3f}".format(ablation["Cohen's kappa"][2])).replace("-", "−")),
        ("Δ Specificity on ablation",  ("{:+.3f}".format(ablation["Specificity"][2])).replace("-", "−")),
    ]
    # Card pinned high so it sits close to the bar chart above.
    _two_col_card(ax_box, pairs,
                  x=0.02, y_top=1.10, width=0.96, line_h=0.135,
                  font_size=8.6, label_frac=0.66)
    # Italic note pulled tight to the BOTTOM EDGE of the card so it reads
    # as the card's own caption, not as Panel C content below.
    ax_box.text(0.03, 0.20,
                "Stable operating point; no-SS ablation reduced specificity.",
                ha="left", va="top", fontsize=7.8, color=GREY_500,
                style="italic", clip_on=False, transform=ax_box.transAxes)


# ──────────────────────────────────────────────────────────────────────────────
# Panel C — Custom HMM bit score vs Tier 2 score (RefSeq n=18 deployment)
# ──────────────────────────────────────────────────────────────────────────────
def panel_c(ax, c_df: pd.DataFrame, cal_summary: dict, ref_summary: dict):
    """Scatter of all 18 calibration-set-novel domains: HMM bit score vs Tier 2 score.

    Strategic labels mark: (i) the two top-bit-score Halorussus EXCLUDED via
    HG1 failure, (ii) the AS7092-supported Q-anchor singleton, and (iii) the
    lowest-bit-score score-positive partial-domain.
    """
    n = len(c_df)

    # Title (short) + subtitle
    ax.set_title("C   HMM score ≠ Tier 2 priority",
                 fontsize=12, fontweight="bold", loc="left",
                 pad=22, color=NAVY_700)
    ax.text(0.0, 1.015,
            f"RefSeq deployment · n = {n} calibration-set-novel domains",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9, style="italic", color=GREY_500)
    # HMM name: "Custom HalC8_v2.2 HMM bit score" — used consistently in
    # the x-axis label below.

    # Y normalized to /7.0 (max Tier 2 score in no-SS mode)
    y_full = c_df["tier2_score"] / 7.0
    c_df = c_df.copy()
    c_df["y"] = y_full

    # Split by score-positive and Halorussus
    sp_other = c_df[(c_df["is_score_positive"] == 1) & (~c_df["is_halorussus"])]
    sp_halo  = c_df[(c_df["is_score_positive"] == 1) & (c_df["is_halorussus"])]
    np_other = c_df[(c_df["is_score_positive"] == 0) & (~c_df["is_halorussus"])]
    np_halo  = c_df[(c_df["is_score_positive"] == 0) & (c_df["is_halorussus"])]

    # Horizontal reference line at score = 3.5 (= y=0.5) — score threshold.
    # zorder=11 places it ABOVE the threshold-note's white bbox (which is
    # at the text's zorder=10) so the dashed line is never cut visually.
    ax.axhline(0.5, color=GREY_500, linestyle="--", linewidth=1.0,
               alpha=0.85, zorder=11)
    # Compact label, anchored at the left edge of the plot, no bbox so the
    # dashed line crosses cleanly through the gaps between characters.
    ax.text(c_df["hmm_bitscore"].min() - 5, 0.49,
            "score ≥ 3.5 + hard gates",
            fontsize=8, color=GREY_500, ha="left", va="top",
            style="italic", zorder=10)

    # Plot points
    if not sp_other.empty:
        ax.scatter(sp_other["hmm_bitscore"], sp_other["y"],
                   s=85, marker="o", facecolor=NAVY_700,
                   edgecolor=WHITE, linewidths=0.9, alpha=0.85, zorder=4,
                   label="Operational score-positive")
    if not np_other.empty:
        ax.scatter(np_other["hmm_bitscore"], np_other["y"],
                   s=85, marker="o", facecolor=WHITE, edgecolor=GREY_500,
                   linewidths=1.5, alpha=0.95, zorder=4,
                   label="Not prioritized")
    if not sp_halo.empty:
        ax.scatter(sp_halo["hmm_bitscore"], sp_halo["y"],
                   s=120, marker="^", facecolor=NAVY_700,
                   edgecolor=WHITE, linewidths=0.9, alpha=0.85, zorder=5)
    if not np_halo.empty:
        ax.scatter(np_halo["hmm_bitscore"], np_halo["y"],
                   s=130, marker="^", facecolor=WHITE, edgecolor=NAVY_900,
                   linewidths=1.8, alpha=0.95, zorder=6,
                   label="Halorussus (HG1 fail)")

    # ── Strategic labels — short, sit in empty quadrants ───────────────
    annotations = []
    indexed = c_df.set_index("accession", verify_integrity=True)
    def point(accession):
        row = indexed.loc[accession]
        return (row["hmm_bitscore"], row["y"])


    # Top-bit-score Halorussus pair — merged label.
    # Genus name italicized via mathtext for journal convention.
    # Placed lower-left (x=70-80, y=0.10-0.26 region) per spec.
    annotations.append(dict(
        text="Top HMM scores  $\\it{Halorussus}$\nEXCLUDED (HG1 fail)",
        xy_data=[point("WP_276279820.1"), point("WP_253478286.1")],
        xy_text=(76, 0.18),
        ha="center", va="center",
    ))

    # AS7092-supported Q-anchor singleton — placed RIGHT of the data cluster
    # at mid-height so it doesn't compete with the inset or with the
    # Halorussus label.
    annotations.append(dict(
        text="WP_015299281.1\nAS7092-supported Q-anchor\nsingleton (Cluster 3)",
        xy_data=[point("WP_015299281.1")],
        xy_text=(92, 0.82),
        ha="center", va="center",
    ))

    # Lowest-bit-score score-positive (partial 50 aa) — LOWER-LEFT, lowered
    # further on Y so it sits just above the x-axis.
    annotations.append(dict(
        text=("WP_010904171.1\nPartial 50-aa (audit only)\n"
              "lowest HMM among score-positive"),
        xy_data=[point("WP_010904171.1")],
        xy_text=(50, 0.04),
        ha="left", va="bottom",
    ))

    for ann in annotations:
        x_text, y_text = ann["xy_text"]
        ax.text(x_text, y_text, ann["text"],
                fontsize=7.6, color=NAVY_900,
                ha=ann["ha"], va=ann["va"],
                bbox=dict(boxstyle="round,pad=0.35", facecolor=NAVY_050,
                          edgecolor=NAVY_300, linewidth=0.7),
                zorder=8)
        for (xd, yd) in ann["xy_data"]:
            ax.annotate("", xy=(xd, yd), xytext=(x_text, y_text),
                        arrowprops=dict(arrowstyle="-", color=NAVY_500,
                                        lw=0.6, alpha=0.7,
                                        connectionstyle="arc3,rad=0.0"),
                        zorder=7)

    # ── Inset metrics box (calibration n=42 + RefSeq n=18) ────────────
    # Placed UPPER-RIGHT, ONE STEP DOWN from the subtitle so it stacks under
    # it without colliding with the WP_010904171.1 label at upper-left.
    # Restored to upper-LEFT (original placement) per user request.
    inset_text = (
        f"Calibration n = {int(cal_summary['n'])}:   AUROC = {float(cal_summary['auroc_hmm_bitscore_vs_reference_binary']):.3f},  best κ = {float(cal_summary['best_hmm_threshold_kappa']):.3f}\n"
        f"RefSeq n = {int(ref_summary['n'])}:         AUROC = {float(ref_summary['auroc_hmm_bitscore_vs_tier2_binary']):.3f},  best κ = {float(ref_summary['best_hmm_threshold_kappa']):.3f},  ρ = {float(ref_summary['spearman_hmm_bitscore_vs_tier2_score']):.3f}"
    )
    ax.text(0.015, 0.97, inset_text, transform=ax.transAxes,
            ha="left", va="top", fontsize=8.0, color=INK,
            family="monospace",
            bbox=dict(boxstyle="round,pad=0.45", facecolor=NAVY_050,
                      edgecolor=NAVY_300, linewidth=0.8),
            zorder=9)

    # (Bottom-right conclusion label removed per spec — title already
    #  carries the message "HMM score ≠ Tier 2 priority".)

    # Axes — same HMM name as subtitle; Y label flags normalization
    ax.set_xlabel("Custom HalC8_v2.2 HMM bit score",
                  fontsize=9.5, fontweight="bold", color=NAVY_700)
    ax.set_ylabel("Tier 2 score / 7.0, no-SS mode",
                  fontsize=9.5, fontweight="bold", color=NAVY_700)
    x_min = c_df["hmm_bitscore"].min()
    x_max = c_df["hmm_bitscore"].max()
    x_pad = (x_max - x_min) * 0.10
    ax.set_xlim(x_min - x_pad - 4, x_max + x_pad + 14)
    ax.set_ylim(0.0, 1.30)
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(True, alpha=0.25, linewidth=0.5, linestyle=":")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # Legend — lower-right, away from the strategic labels
    ax.legend(loc="lower right", fontsize=8.0, frameon=True,
              framealpha=0.95, edgecolor=NAVY_300, borderpad=0.5,
              handletextpad=0.5)


# ──────────────────────────────────────────────────────────────────────────────
# Composer
# ──────────────────────────────────────────────────────────────────────────────
def build_figure(a, b_summary, ablation, c_df, cal_summary, ref_summary) -> plt.Figure:
    """Layout: A (top-left) + B (top-right) + C (full-width bottom).

    Sized at BMC Bioinformatics print dimensions (170 mm × 200 mm) so the
    figure renders at its publication size without downscaling — Panel C
    in particular stays legible after editorial typesetting.
    """
    fig = plt.figure(figsize=(170 / 25.4, 200 / 25.4), facecolor="white")

    # Three rows × two cols. Top row = matrix + chart, mid row = cards, bot row = C.
    gs = GridSpec(
        nrows=3, ncols=2, figure=fig,
        width_ratios=[1.0, 1.0],
        height_ratios=[1.4, 0.82, 1.7],
        wspace=0.30, hspace=0.35,
        left=0.085, right=0.975, top=0.955, bottom=0.055,
    )

    ax_a_plot = fig.add_subplot(gs[0, 0])
    ax_a_box  = fig.add_subplot(gs[1, 0])
    ax_b_plot = fig.add_subplot(gs[0, 1])
    ax_b_box  = fig.add_subplot(gs[1, 1])
    ax_c      = fig.add_subplot(gs[2, :])

    panel_a(ax_a_plot, ax_a_box, a)
    panel_b(ax_b_plot, ax_b_box, b_summary, ablation)
    panel_c(ax_c, c_df, cal_summary, ref_summary)

    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
def main():
    here = Path(__file__).resolve().parent
    default_indir  = here / "inputs"
    default_outdir = here

    parser = argparse.ArgumentParser(description="Generate Figure 2 from frozen TSV inputs.")
    parser.add_argument("--indir",  default=str(default_indir),
                        help="Directory containing the frozen Figure 2 TSV inputs.")
    parser.add_argument("--outdir", default=str(default_outdir),
                        help="Output directory.")
    args = parser.parse_args()

    indir  = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("generate_figure2.py — HalC8 manuscript Figure 2")
    print("=" * 72)
    print(f"  Input dir : {indir}")
    print(f"  Output dir: {outdir}")
    print()

    a         = load_panel_a       (indir / "figure2_panelA_metrics.tsv")
    b_summary = load_panel_b_summary(indir / "figure2_panelB_summary.tsv")
    ablation = load_panel_b_ablation(indir / "figure2_panelB_no_ss_ablation.tsv")
    cal_summary = _read_wide_tsv(indir / "calibration_42_summary.tsv")
    ref_summary = _read_wide_tsv(indir / "refseq_18_summary.tsv")
    c_df      = load_panel_c       (indir / "figure2_panelC_hmm_vs_tier2.tsv")

    print(f"  Panel A: n={a['n']} TP/TN/FP/FN={a['tp']}/{a['tn']}/{a['fp']}/{a['fn']}  κ={a['kappa']:.3f}")
    print(f"  Panel B: fixed operating point (ML={b_summary['ml']}, PROB={b_summary['pr']}) "
          f"{b_summary['n_agree']}/{b_summary['n_folds']} folds  κ={b_summary['kappa']:.3f}")
    print(f"  Panel B: calibration-frame ablation (Table S9 Panel C)")
    print(f"  Panel C: {len(c_df)} RefSeq domains")
    if len(c_df) != int(ref_summary['n']) or a['n'] != int(cal_summary['n']):
        raise ValueError("Input sample counts disagree")
    for key, name in [("Cohen's kappa", "kappa"), ("Sensitivity", "sensitivity"), ("Specificity", "specificity"), ("Precision", "precision"), ("F1", "f1")]:
        if abs(ablation[key][0] - a[name]) > 0.0005:
            raise ValueError(f"Panel A and B disagree: {key}")
    fig = build_figure(a, b_summary, ablation, c_df, cal_summary, ref_summary)

    # Final manuscript ordering:
    #   Figure 1 = workflow schematic
    #   Figure 2 = calibration / no-SS / HMM-only comparator   ← this script
    #   Figure 3 = RefSeq deployment funnel
    #   Figure 4 = cluster-level synteny / triage
    svg_path = outdir / "Figure_2.svg"
    pdf_path = outdir / "Figure_2.pdf"
    fig.savefig(svg_path, format="svg", bbox_inches="tight",
                facecolor="white", edgecolor="none")
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight",
                facecolor="white", edgecolor="none")
    fig.savefig(outdir / "Figure_2.png", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"✓ SVG saved: {svg_path}")
    print(f"✓ PDF saved: {pdf_path}")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
