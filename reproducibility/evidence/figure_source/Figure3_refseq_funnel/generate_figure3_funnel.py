#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_figure3_funnel.py
═══════════════════════════════════════════════════════════════════════════════
Figure 3 — HalC8-family candidate-domain prioritization funnel
(Halobacteria RefSeq → frozen v2.2 binary scoring → cluster de-duplication)

Final manuscript figure ordering:
  Figure 1 = workflow schematic
  Figure 2 = calibration / no-SS / HMM-only comparator
  Figure 3 = RefSeq deployment funnel        ← this script
  Figure 4 = cluster-level synteny / triage

DESIGN RULES (locked):
  • BINARY framing.  The visual narrative is 18 calibration-set-novel domains →
    (13 operational score-positive outputs + 5 not-prioritized domains), with
    the score-positive branch explicitly split into 12 fully assessable domains
    plus one Z3-imputed below-calibration-range partial-domain output. Cluster
    de-duplication is interpreted as 5 operational clusters, comprising 4 fully
    assessable clusters plus Cluster 5 retained for audit/provenance only.
  • The legacy 4-tier MOST LIKELY / PROBABLE / MARGINAL / EXCLUDED
    hierarchy is NEVER visualised as primary bars or panels — it
    appears only as small "audit labels: …" annotations inside each
    binary box, and in the footer * mapping, for reproducibility
    traceability only.
  • No figure title inside the image (BMC requires titles in the
    figure legend, not on the artwork itself).
  • All numeric values + label strings are read from a frozen TSV.
  • Output is a single editable SVG plus a vectorial PDF (Helvetica,
    8–12 pt body type, 170 mm canvas).

INPUT  (default: ./inputs/Figure3.tsv when run from Figure3_refseq_funnel_source/)
OUTPUTS  (default: ./figures/)
  Figure_3.svg   — editable, live text, vector icons
  Figure_3.pdf   — vector PDF (cairosvg, rsvg-convert, or inkscape)

USAGE
  python3 generate_figure3_funnel.py
  python3 generate_figure3_funnel.py --in <tsv> --out <dir>
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
from xml.sax.saxutils import escape


# ──────────────────────────────────────────────────────────────────────────────
# Palette  (kept close to the reference target PDF, BMC-friendly contrasts)
# ──────────────────────────────────────────────────────────────────────────────
COL_NAVY        = "#1f3a5f"
COL_NAVY_FILL   = "#dfeaf5"   # tier 1 trapezoid fill
COL_NAVY_LINE   = "#6f8fb5"   # tier 1 trapezoid stroke
COL_NAVY_BIG    = "#1a3358"
COL_NAVY_SUB    = "#6f8fb5"

COL_TEAL        = "#1f7a5e"
COL_TEAL_FILL   = "#dff0e8"   # tier 2 trapezoid fill
COL_TEAL_LINE   = "#6fa890"
COL_TEAL_BIG    = "#1a5c47"
COL_TEAL_SUB    = "#5a8a78"

COL_GREEN       = "#2c8a3a"   # tier 3 + score-positive box
COL_GREEN_FILL  = "#e6f3da"
COL_GREEN_LINE  = "#88b870"
COL_GREEN_BIG   = "#1a5c2a"
COL_GREEN_TXT   = "#2a6a38"
COL_GREEN_AUDIT = "#7aa088"   # very subtle, audit-label gray-green

COL_PEACH_FILL  = "#fce0c8"   # not-prioritized box
COL_PEACH_LINE  = "#d49060"
COL_PEACH_BIG   = "#a04018"
COL_PEACH_TXT   = "#9a4828"
COL_PEACH_AUDIT = "#b88060"

COL_CLUSTER_FILL = "#d8efe6"  # cluster-level de-duplication band
COL_CLUSTER_LINE = "#5a9a88"
COL_CLUSTER_TXT  = "#1a5c4a"

COL_INPUT_FILL  = "#eef4fb"   # input proteome box
COL_INPUT_LINE  = "#a4bcd8"

COL_REDUCTION_LINE = "#1f3a5f"
COL_REDUCTION_TXT  = "#1f3a5f"

COL_AXIS_GREY   = "#6e6e6e"
COL_FOOTER_TXT  = "#444444"

WHITE = "#ffffff"


# ──────────────────────────────────────────────────────────────────────────────
# I/O helpers
# ──────────────────────────────────────────────────────────────────────────────
def read_tsv(path: Path) -> dict:
    if not path.is_file():
        sys.exit(f"ERROR: input TSV not found: {path}")
    d = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.rstrip("\n")
            if not s.strip() or s.startswith("#") or "\t" not in s:
                continue
            k, v = s.split("\t", 1)
            k = k.strip().lower()
            v = v.strip()
            if k in ("key", "metric", "field", "name", "parameter"):
                continue
            d[k] = v
    return d


def _i(d: dict, k: str) -> int:
    if k not in d:
        sys.exit(f"ERROR: missing key in TSV: {k}")
    try:
        return int(round(float(d[k])))
    except Exception:
        sys.exit(f"ERROR: TSV key '{k}' is not numeric: {d[k]!r}")


def _f(d: dict, k: str) -> float:
    if k not in d:
        sys.exit(f"ERROR: missing key in TSV: {k}")
    try:
        return float(d[k])
    except Exception:
        sys.exit(f"ERROR: TSV key '{k}' is not numeric: {d[k]!r}")


def _s(d: dict, k: str) -> str:
    if k not in d:
        sys.exit(f"ERROR: missing key in TSV: {k}")
    return d[k]


def _comma(n: int) -> str:
    return f"{n:,}"


def e(s: str) -> str:
    return escape(s, {'"': "&quot;"})


# ──────────────────────────────────────────────────────────────────────────────
# Icon helpers — small monoline SVG glyphs, centred at (0, 0) of a <g transform>
# ──────────────────────────────────────────────────────────────────────────────
def icon_database(stroke: str) -> str:
    """Small DB cylinder, ~18 px tall."""
    return (
        f'<ellipse cx="0" cy="-7" rx="8" ry="2.5" fill="none" stroke="{stroke}" stroke-width="1.4"/>'
        f'<path d="M -8 -7 L -8 7 A 8 2.5 0 0 0 8 7 L 8 -7" fill="none" stroke="{stroke}" stroke-width="1.4"/>'
        f'<path d="M -8 -1 A 8 2.5 0 0 0 8 -1" fill="none" stroke="{stroke}" stroke-width="1.4"/>'
        f'<path d="M -8 5 A 8 2.5 0 0 0 8 5" fill="none" stroke="{stroke}" stroke-width="1.4"/>'
    )


def icon_lookup(stroke: str) -> str:
    """Magnifying glass."""
    return (
        f'<circle cx="-2" cy="-2" r="5.5" fill="none" stroke="{stroke}" stroke-width="1.5"/>'
        f'<line x1="2.4" y1="2.4" x2="7.5" y2="7.5" stroke="{stroke}" stroke-width="1.7" stroke-linecap="round"/>'
    )


def icon_clipboard(stroke: str) -> str:
    """Audit clipboard with two check marks."""
    return (
        f'<rect x="-6" y="-7" width="12" height="15" rx="1.4" fill="none" stroke="{stroke}" stroke-width="1.4"/>'
        f'<rect x="-3" y="-9" width="6" height="3" rx="0.6" fill="none" stroke="{stroke}" stroke-width="1.4"/>'
        f'<path d="M -4 -2 L -2 0 L 2 -3" fill="none" stroke="{stroke}" stroke-width="1.3" stroke-linecap="round" stroke-linejoin="round"/>'
        f'<path d="M -4 3 L -2 5 L 2 2" fill="none" stroke="{stroke}" stroke-width="1.3" stroke-linecap="round" stroke-linejoin="round"/>'
    )


def icon_dial(stroke: str) -> str:
    """Scoring dial / speedometer."""
    return (
        f'<path d="M -8 4 A 8 8 0 0 1 8 4" fill="none" stroke="{stroke}" stroke-width="1.5"/>'
        f'<line x1="0" y1="4" x2="4" y2="-3" stroke="{stroke}" stroke-width="1.7" stroke-linecap="round"/>'
        f'<circle cx="0" cy="4" r="1.4" fill="{stroke}"/>'
    )


def icon_shield_check(stroke: str) -> str:
    """Shield with tick (score-positive)."""
    return (
        f'<path d="M -7 -7 L 7 -7 L 7 1 Q 7 7 0 9 Q -7 7 -7 1 Z" '
        f'fill="none" stroke="{stroke}" stroke-width="1.5" stroke-linejoin="round"/>'
        f'<path d="M -3 -1 L -0.5 2 L 3.5 -3" fill="none" '
        f'stroke="{stroke}" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/>'
    )


def icon_shield_x(stroke: str) -> str:
    """Shield with × (not prioritized)."""
    return (
        f'<path d="M -7 -7 L 7 -7 L 7 1 Q 7 7 0 9 Q -7 7 -7 1 Z" '
        f'fill="none" stroke="{stroke}" stroke-width="1.5" stroke-linejoin="round"/>'
        f'<line x1="-2.5" y1="-3" x2="2.5" y2="2" stroke="{stroke}" stroke-width="1.6" stroke-linecap="round"/>'
        f'<line x1="-2.5" y1="2"  x2="2.5" y2="-3" stroke="{stroke}" stroke-width="1.6" stroke-linecap="round"/>'
    )


def icon_cluster(stroke: str) -> str:
    """Three interconnected nodes — cluster de-duplication."""
    return (
        f'<line x1="-5" y1="-3" x2="5" y2="-3" stroke="{stroke}" stroke-width="1.3"/>'
        f'<line x1="-5" y1="-3" x2="0" y2="6"  stroke="{stroke}" stroke-width="1.3"/>'
        f'<line x1="5"  y1="-3" x2="0" y2="6"  stroke="{stroke}" stroke-width="1.3"/>'
        f'<circle cx="-5" cy="-3" r="2.6" fill="{WHITE}" stroke="{stroke}" stroke-width="1.5"/>'
        f'<circle cx="5"  cy="-3" r="2.6" fill="{WHITE}" stroke="{stroke}" stroke-width="1.5"/>'
        f'<circle cx="0"  cy="6"  r="2.6" fill="{WHITE}" stroke="{stroke}" stroke-width="1.5"/>'
    )


# ──────────────────────────────────────────────────────────────────────────────
# Numbered marker — filled circle with white digit, used at upper-left of each
# tier trapezoid as a visual stamp.
# ──────────────────────────────────────────────────────────────────────────────
def numbered_marker(cx: int, cy: int, n: int, fill: str, r: int = 14) -> str:
    return (
        f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}" stroke="{WHITE}" stroke-width="1.6"/>'
        f'<text x="{cx}" y="{cy + 4.5}" text-anchor="middle" fill="{WHITE}" '
        f'font-size="14" font-weight="bold">{n}</text>'
    )


# ──────────────────────────────────────────────────────────────────────────────
# Trapezoid path helper.  Given top and bottom edge x-ranges plus y bounds,
# returns the SVG `d` attribute for a slightly rounded trapezoid.
# ──────────────────────────────────────────────────────────────────────────────
def trapezoid(x_top_l: int, x_top_r: int,
              x_bot_l: int, x_bot_r: int,
              y_top: int, y_bot: int, r: int = 4) -> str:
    # Plain trapezoid with small rounded corners on the top edge only.
    return (
        f"M {x_top_l + r},{y_top} "
        f"L {x_top_r - r},{y_top} "
        f"Q {x_top_r},{y_top} {x_top_r},{y_top + r} "
        f"L {x_bot_r},{y_bot - r} "
        f"Q {x_bot_r},{y_bot} {x_bot_r - r},{y_bot} "
        f"L {x_bot_l + r},{y_bot} "
        f"Q {x_bot_l},{y_bot} {x_bot_l},{y_bot - r} "
        f"L {x_top_l},{y_top + r} "
        f"Q {x_top_l},{y_top} {x_top_l + r},{y_top} Z"
    )


# ──────────────────────────────────────────────────────────────────────────────
# SVG builder
# ──────────────────────────────────────────────────────────────────────────────
def build_svg(d: dict) -> str:

    # ── Numbers + labels from TSV ────────────────────────────────────────
    proteome_n         = _i(d, "proteome_n")
    proteome_label     = _s(d, "proteome_label")
    hmm_positive_n     = _i(d, "hmm_positive_n")
    calibration_set_novel_n   = _i(d, "calibration_set_novel_n")
    score_positive_n   = _i(d, "score_positive_n")
    not_prioritized_n  = _i(d, "not_prioritized_n")
    fully_assessable_scorepositive_n = _i(d, "fully_assessable_scorepositive_n")
    z3_imputed_partial_scorepositive_n = _i(d, "z3_imputed_partial_scorepositive_n")
    fully_assessable_cluster_n = _i(d, "fully_assessable_cluster_n")
    z3_imputed_partial_cluster_n = _i(d, "z3_imputed_partial_cluster_n")


    tier1_label        = _s(d, "tier1_label")
    tier1_subtitle     = _s(d, "tier1_subtitle")
    novelty_audit_label = _s(d, "novelty_audit_label")
    novelty_audit_subtitle = _s(d, "novelty_audit_subtitle")
    tier2_label        = _s(d, "tier2_label")
    tier2_subtitle     = _s(d, "tier2_subtitle")

    score_positive_subtitle  = _s(d, "score_positive_subtitle")
    not_prioritized_subtitle = _s(d, "not_prioritized_subtitle")
    audit_label_score_positive  = _s(d, "audit_label_score_positive")
    audit_label_not_prioritized = _s(d, "audit_label_not_prioritized")

    cluster_block_title = _s(d, "cluster_block_title")
    cluster_block_body  = _s(d, "cluster_block_body")

    footer_main    = _s(d, "footer_main")
    footer_partial = _s(d, "footer_partial")

    # Fail closed if frozen flow counts do not reconcile.
    ml, pr, mg, ex = (_i(d, k) for k in ("most_likely_n", "probable_n", "marginal_n", "excluded_n"))
    if not (hmm_positive_n - _i(d, "calibration_removed_n") == calibration_set_novel_n
            and ml + pr == score_positive_n and mg + ex == not_prioritized_n
            and score_positive_n + not_prioritized_n == calibration_set_novel_n
            and fully_assessable_scorepositive_n + z3_imputed_partial_scorepositive_n == score_positive_n
            and fully_assessable_cluster_n + z3_imputed_partial_cluster_n == _i(d, "clusters_n")):
        raise ValueError("Frozen Figure 3 counts do not reconcile")

    # ── Canvas ──────────────────────────────────────────────────────────
    # The viewBox is the internal coordinate grid (790×720 user units).
    # Declaring an explicit physical width="170mm" pins the rendered PDF
    # to BMC Bioinformatics' double-column width so that every SVG
    # font-size value translates literally to printed points (1 SVG
    # unit ≈ 1 pt at that physical size). Without this attribute,
    # cairosvg defaults to 96 DPI and produces a 209 mm-wide PDF that
    # BMC's typesetter would scale to 170 mm, silently shrinking all
    # body text to ~5–7 pt and breaking the 8 pt minimum.
    W, H = 790, 720
    BMC_WIDTH_MM  = 170.0
    BMC_HEIGHT_MM = BMC_WIDTH_MM * H / W  # preserve aspect ratio

    parts: list[str] = []
    parts.append(
        f'<svg viewBox="0 0 {W} {H}" '
        f'width="{BMC_WIDTH_MM:.2f}mm" height="{BMC_HEIGHT_MM:.2f}mm" '
        f'xmlns="http://www.w3.org/2000/svg" '
        f'role="img" font-family="Helvetica, Arial, sans-serif">'
    )
    parts.append(
        f'  <title>HalC8-family candidate-domain prioritization funnel</title>'
    )
    parts.append(f'  <rect width="{W}" height="{H}" fill="{WHITE}"/>')

    # ─────────────────────────────────────────────────────────────────
    # 1. INPUT PROTEOME — wide rounded box, navy outline
    # ─────────────────────────────────────────────────────────────────
    parts.append(
        f'  <rect x="60" y="20" width="670" height="80" rx="10" '
        f'fill="{COL_INPUT_FILL}" stroke="{COL_NAVY_LINE}" stroke-width="1.2"/>'
    )
    parts.append(
        f'  <g transform="translate(100,60)">{icon_database(COL_NAVY)}</g>'
    )
    parts.append(
        f'  <text x="395" y="44" text-anchor="middle" fill="{COL_NAVY}" '
        f'font-size="11.5" font-weight="bold">Input proteome set</text>'
    )
    parts.append(
        f'  <text x="395" y="75" text-anchor="middle" fill="{COL_NAVY_BIG}" '
        f'font-size="26" font-weight="bold">{_comma(proteome_n)}</text>'
    )
    parts.append(
        f'  <text x="395" y="93" text-anchor="middle" fill="{COL_NAVY_SUB}" '
        f'font-size="9.5" font-style="italic">{e(proteome_label)}</text>'
    )

    # ─────────────────────────────────────────────────────────────────
    # FUNNEL — three tier trapezoids, each narrower than the previous
    # ─────────────────────────────────────────────────────────────────
    # Layout coordinates (chosen so successive trapezoids share an edge)
    t1_top_l, t1_top_r = 60,  730
    t1_bot_l, t1_bot_r = 130, 660
    t1_y_top, t1_y_bot = 115, 220

    t2_top_l, t2_top_r = t1_bot_l, t1_bot_r
    t2_bot_l, t2_bot_r = 200, 590
    t2_y_top, t2_y_bot = 230, 330

    t3_top_l, t3_top_r = t2_bot_l, t2_bot_r
    t3_bot_l, t3_bot_r = 270, 520
    t3_y_top, t3_y_bot = 340, 445

    # ── Tier 1 trapezoid (navy fill) ──────────────────────────────────
    parts.append(
        f'  <path d="{trapezoid(t1_top_l, t1_top_r, t1_bot_l, t1_bot_r, t1_y_top, t1_y_bot)}" '
        f'fill="{COL_NAVY_FILL}" stroke="{COL_NAVY_LINE}" stroke-width="1.2"/>'
    )
    parts.append(numbered_marker(95, 130, 1, COL_NAVY))
    parts.append(f'  <g transform="translate(135,130)">{icon_lookup(COL_NAVY)}</g>')
    parts.append(
        f'  <text x="400" y="135" text-anchor="middle" fill="{COL_NAVY}" '
        f'font-size="11" font-weight="bold">{e(tier1_label)}</text>'
    )
    parts.append(
        f'  <text x="400" y="151" text-anchor="middle" fill="{COL_NAVY_SUB}" '
        f'font-size="9" font-style="italic">{e(tier1_subtitle)}</text>'
    )
    parts.append(
        f'  <text x="400" y="195" text-anchor="middle" fill="{COL_NAVY_BIG}" '
        f'font-size="28" font-weight="bold">{hmm_positive_n}</text>'
    )
    parts.append(
        f'  <text x="400" y="214" text-anchor="middle" fill="{COL_NAVY_BIG}" '
        f'font-size="10.5">HMM-positive domains</text>'
    )

    # ── Tier 2 trapezoid (teal fill) ──────────────────────────────────
    parts.append(
        f'  <path d="{trapezoid(t2_top_l, t2_top_r, t2_bot_l, t2_bot_r, t2_y_top, t2_y_bot)}" '
        f'fill="{COL_TEAL_FILL}" stroke="{COL_TEAL_LINE}" stroke-width="1.2"/>'
    )
    parts.append(numbered_marker(150, 245, 2, COL_TEAL))
    parts.append(f'  <g transform="translate(190,245)">{icon_clipboard(COL_TEAL)}</g>')
    parts.append(
        f'  <text x="400" y="250" text-anchor="middle" fill="{COL_TEAL_BIG}" '
        f'font-size="11" font-weight="bold">{e(novelty_audit_label)}</text>'
    )
    parts.append(
        f'  <text x="400" y="266" text-anchor="middle" fill="{COL_TEAL_SUB}" '
        f'font-size="9" font-style="italic">{e(novelty_audit_subtitle)}</text>'
    )
    parts.append(
        f'  <text x="400" y="305" text-anchor="middle" fill="{COL_TEAL_BIG}" '
        f'font-size="26" font-weight="bold">{calibration_set_novel_n}</text>'
    )
    parts.append(
        f'  <text x="400" y="323" text-anchor="middle" fill="{COL_TEAL_BIG}" '
        f'font-size="10.5">calibration-set-novel domains</text>'
    )

    # ── Tier 3 trapezoid (green fill) ─────────────────────────────────
    parts.append(
        f'  <path d="{trapezoid(t3_top_l, t3_top_r, t3_bot_l, t3_bot_r, t3_y_top, t3_y_bot)}" '
        f'fill="{COL_GREEN_FILL}" stroke="{COL_GREEN_LINE}" stroke-width="1.2"/>'
    )
    parts.append(numbered_marker(220, 355, 3, COL_GREEN))
    parts.append(f'  <g transform="translate(260,355)">{icon_dial(COL_GREEN)}</g>')
    parts.append(
        f'  <text x="400" y="360" text-anchor="middle" fill="{COL_GREEN_BIG}" '
        f'font-size="11" font-weight="bold">{e(tier2_label)}</text>'
    )
    parts.append(
        f'  <text x="400" y="376" text-anchor="middle" fill="{COL_GREEN_TXT}" '
        f'font-size="9" font-style="italic">{e(tier2_subtitle)}</text>'
    )
    parts.append(
        f'  <text x="400" y="418" text-anchor="middle" fill="{COL_GREEN_BIG}" '
        f'font-size="26" font-weight="bold">{calibration_set_novel_n}</text>'
    )
    parts.append(
        f'  <text x="400" y="436" text-anchor="middle" fill="{COL_GREEN_BIG}" '
        f'font-size="10.5">domains scored</text>'
    )

    # ─────────────────────────────────────────────────────────────────
    # 5. BIFURCATION arrows — solid green → score-positive,
    #                          dashed grey → not-prioritized
    # ─────────────────────────────────────────────────────────────────
    # Solid arrow: from bottom of tier 3 down-left to score-positive box
    parts.append(
        f'  <path d="M 365 448 L 365 470 L 230 470 L 230 488" '
        f'fill="none" stroke="{COL_GREEN}" stroke-width="1.8" '
        f'stroke-linecap="round" stroke-linejoin="round"/>'
    )
    parts.append(
        f'  <polygon points="230,498 225,488 235,488" fill="{COL_GREEN}"/>'
    )
    # Dashed arrow: from bottom of tier 3 down-right to not-prioritized box
    parts.append(
        f'  <path d="M 425 448 L 425 470 L 555 470 L 555 488" '
        f'fill="none" stroke="{COL_PEACH_LINE}" stroke-width="1.5" '
        f'stroke-linecap="round" stroke-linejoin="round" stroke-dasharray="5,4"/>'
    )
    parts.append(
        f'  <polygon points="555,498 550,488 560,488" fill="{COL_PEACH_LINE}"/>'
    )

    # ─────────────────────────────────────────────────────────────────
    # 6. TWO PARALLEL BOXES — the visual centre of the figure
    # ─────────────────────────────────────────────────────────────────
    # Score-positive (green) — LEFT
    sp_x, sp_y, sp_w, sp_h = 70, 495, 340, 110
    parts.append(
        f'  <rect x="{sp_x}" y="{sp_y}" width="{sp_w}" height="{sp_h}" rx="8" '
        f'fill="{COL_GREEN_FILL}" stroke="{COL_GREEN_LINE}" stroke-width="1.4"/>'
    )
    parts.append(
        f'  <g transform="translate({sp_x + 28}, {sp_y + 38})">{icon_shield_check(COL_GREEN)}</g>'
    )
    parts.append(
        f'  <text x="{sp_x + 50}" y="{sp_y + 34}" text-anchor="start" fill="{COL_GREEN_BIG}" '
        f'font-size="12.5" font-weight="bold">{score_positive_n} operational score-positive outputs</text>'
    )
    parts.append(
        f'  <text x="{sp_x + 50}" y="{sp_y + 51}" text-anchor="start" fill="{COL_GREEN_TXT}" '
        f'font-size="9.5">{e(score_positive_subtitle)}</text>'
    )
    parts.append(
        f'  <rect x="{sp_x + 50}" y="{sp_y + 61}" width="145" height="21" rx="5" '
        f'fill="{WHITE}" stroke="{COL_GREEN_LINE}" stroke-width="0.8" opacity="0.95"/>'
    )
    parts.append(
        f'  <text x="{sp_x + 122.5}" y="{sp_y + 75}" text-anchor="middle" fill="{COL_GREEN_BIG}" '
        f'font-size="8.6" font-weight="bold">{fully_assessable_scorepositive_n} fully assessable</text>'
    )
    parts.append(
        f'  <rect x="{sp_x + 205}" y="{sp_y + 61}" width="118" height="21" rx="5" '
        f'fill="{WHITE}" stroke="{COL_GREEN_LINE}" stroke-width="0.8" stroke-dasharray="3,3" opacity="0.95"/>'
    )
    parts.append(
        f'  <text x="{sp_x + 264}" y="{sp_y + 75}" text-anchor="middle" fill="{COL_GREEN_TXT}" '
        f'font-size="8.4" font-weight="bold">{z3_imputed_partial_scorepositive_n} partial*</text>'
    )
    parts.append(
        f'  <line x1="{sp_x + 16}" y1="{sp_y + 86}" x2="{sp_x + sp_w - 16}" y2="{sp_y + 86}" '
        f'stroke="{COL_GREEN_LINE}" stroke-width="0.8" stroke-dasharray="3,3" opacity="0.6"/>'
    )
    parts.append(
        f'  <text x="{sp_x + sp_w/2}" y="{sp_y + 103}" text-anchor="middle" fill="{COL_GREEN_AUDIT}" '
        f'font-size="8.5" font-style="italic">{e(audit_label_score_positive)}</text>'
    )

    # Not-prioritized (peach) — RIGHT
    np_x, np_y, np_w, np_h = 420, 495, 305, 110
    parts.append(
        f'  <rect x="{np_x}" y="{np_y}" width="{np_w}" height="{np_h}" rx="8" '
        f'fill="{COL_PEACH_FILL}" stroke="{COL_PEACH_LINE}" stroke-width="1.4"/>'
    )
    parts.append(
        f'  <g transform="translate({np_x + 28}, {np_y + 38})">{icon_shield_x(COL_PEACH_BIG)}</g>'
    )
    parts.append(
        f'  <text x="{np_x + 50}" y="{np_y + 34}" text-anchor="start" fill="{COL_PEACH_BIG}" '
        f'font-size="14" font-weight="bold">{not_prioritized_n} not-prioritized domains</text>'
    )
    parts.append(
        f'  <text x="{np_x + 50}" y="{np_y + 51}" text-anchor="start" fill="{COL_PEACH_TXT}" '
        f'font-size="9.5">{e(not_prioritized_subtitle)}</text>'
    )
    parts.append(
        f'  <line x1="{np_x + 16}" y1="{np_y + 73}" x2="{np_x + np_w - 16}" y2="{np_y + 73}" '
        f'stroke="{COL_PEACH_LINE}" stroke-width="0.8" stroke-dasharray="3,3" opacity="0.6"/>'
    )
    parts.append(
        f'  <text x="{np_x + np_w/2}" y="{np_y + 92}" text-anchor="middle" fill="{COL_PEACH_AUDIT}" '
        f'font-size="8.5" font-style="italic">{e(audit_label_not_prioritized)}</text>'
    )

    # Vertical arrow from score-positive box down into cluster block
    parts.append(
        f'  <line x1="{sp_x + sp_w/2}" y1="{sp_y + sp_h + 2}" x2="{sp_x + sp_w/2}" y2="{sp_y + sp_h + 16}" '
        f'stroke="{COL_GREEN}" stroke-width="1.6"/>'
    )
    parts.append(
        f'  <polygon points="{sp_x + sp_w/2 - 5},{sp_y + sp_h + 14} '
        f'{sp_x + sp_w/2 + 5},{sp_y + sp_h + 14} '
        f'{sp_x + sp_w/2},{sp_y + sp_h + 22}" fill="{COL_GREEN}"/>'
    )

    # ─────────────────────────────────────────────────────────────────
    # 7. CLUSTER-LEVEL DE-DUPLICATION BLOCK — full-width teal card
    # ─────────────────────────────────────────────────────────────────
    cb_x, cb_y, cb_w, cb_h = 70, 625, 655, 50
    parts.append(
        f'  <rect x="{cb_x}" y="{cb_y}" width="{cb_w}" height="{cb_h}" rx="8" '
        f'fill="{COL_CLUSTER_FILL}" stroke="{COL_CLUSTER_LINE}" stroke-width="1.3"/>'
    )
    parts.append(
        f'  <g transform="translate({cb_x + 32}, {cb_y + cb_h/2})">{icon_cluster(COL_CLUSTER_LINE)}</g>'
    )
    parts.append(
        f'  <text x="{cb_x + 60}" y="{cb_y + 21}" text-anchor="start" fill="{COL_CLUSTER_TXT}" '
        f'font-size="11.5" font-weight="bold">{e(cluster_block_title)}</text>'
    )
    parts.append(
        f'  <text x="{cb_x + 60}" y="{cb_y + 38}" text-anchor="start" fill="{COL_CLUSTER_TXT}" '
        f'font-size="9.5">{fully_assessable_cluster_n} fully assessable clusters + Cluster 5 partial singleton (audit only); clustering threshold {_f(d, "cluster_identity"):.2f}</text>'
    )

    # ─────────────────────────────────────────────────────────────────
    # 8. FOOTERS — two italic notes (* mapping, † partial domain)
    # ─────────────────────────────────────────────────────────────────
    # Footers at 9 pt italic (within BMC 8–12 pt band, a touch more legible
    # at print than the prior 8.5 pt). The mapping `score-positive =
    # ML+PROB / not prioritized = MARG+EXCL` and the audit-traceability
    # disclaimer must remain readable — this is the only place where the
    # legacy 4-tier names appear in canonical form.
    parts.append(
        f'  <text x="60" y="693" text-anchor="start" fill="{COL_FOOTER_TXT}" '
        f'font-size="9" font-style="italic">{e(footer_main)}</text>'
    )
    parts.append(
        f'  <text x="60" y="708" text-anchor="start" fill="{COL_FOOTER_TXT}" '
        f'font-size="9" font-style="italic">{e(footer_partial)}</text>'
    )

    parts.append('</svg>')
    return "\n".join(parts) + "\n"


# ──────────────────────────────────────────────────────────────────────────────
# SVG → PDF
# ──────────────────────────────────────────────────────────────────────────────
def svg_to_pdf(svg_path: Path, pdf_path: Path) -> bool:
    try:
        import cairosvg  # type: ignore
        cairosvg.svg2pdf(url=str(svg_path), write_to=str(pdf_path))
        return True
    except Exception:
        pass
    import shutil, subprocess
    for cmd in (["rsvg-convert", "-f", "pdf", "-o", str(pdf_path), str(svg_path)],
                ["inkscape", str(svg_path), "--export-type=pdf",
                 f"--export-filename={pdf_path}"]):
        if shutil.which(cmd[0]):
            try:
                subprocess.run(cmd, check=True,
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)
                return True
            except Exception:
                continue
    return False


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
def main() -> int:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(
        description="Generate Figure 3 (RefSeq deployment funnel, binary) from frozen TSV."
    )
    p.add_argument("--in",  dest="inp",
                   default=str(here / "inputs" / "Figure3.tsv"))
    p.add_argument("--out", dest="out", default=str(here))
    args = p.parse_args()

    in_path = Path(args.inp).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("generate_figure3_funnel.py — HalC8 main-text Figure 3 (binary funnel)")
    print("=" * 70)
    print(f"  Input : {in_path}")
    print(f"  Output: {out_dir}")

    d = read_tsv(in_path)
    svg = build_svg(d)

    svg_path = out_dir / "Figure_3.svg"
    pdf_path = out_dir / "Figure_3.pdf"
    svg_path.write_text(svg, encoding="utf-8")
    print(f"✓ SVG saved: {svg_path}")

    if svg_to_pdf(svg_path, pdf_path):
        print(f"✓ PDF saved: {pdf_path}")
    else:
        print("ℹ  PDF skipped — install cairosvg (`pip install cairosvg`) "
              "or rsvg-convert / inkscape.")
    import cairosvg
    png_path = out_dir / "Figure_3.png"
    cairosvg.svg2png(url=str(svg_path), write_to=str(png_path), dpi=600)
    from PIL import Image
    with Image.open(png_path) as im:
        im.load()
        im.save(png_path, dpi=(600, 600))
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
