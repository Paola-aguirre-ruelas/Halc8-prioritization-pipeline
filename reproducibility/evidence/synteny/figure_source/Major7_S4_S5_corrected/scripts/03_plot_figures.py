#!/usr/bin/env python3
"""Draw Supplementary Figures S4 and S5 from the visualization annotation tables.

Geometry comes only from annotations/Figure_<panel>_cds_annotation.csv (GenBank coordinates) and
annotations/Figure_<panel>_links.tsv (archived Clinker links). Display normalisation: loci are aligned
at the first base of the candidate CDS (reference: proC8) and loci whose candidate lies on the minus
strand of the archived window are drawn reverse-complemented (marked with a dagger)."""
import csv, json, sys
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, PathPatch, Rectangle
from matplotlib.path import Path as MPath
from matplotlib.lines import Line2D

PKG = Path(__file__).resolve().parent.parent
CFG = json.loads((PKG / "config" / "figure_config.json").read_text())
PAL = {r["display_class"]: r for r in csv.DictReader(open(PKG / "config" / "colour_map.csv"))}
plt.rcParams.update({"font.family": "Arial", "svg.fonttype": "none", "pdf.fonttype": 42,
                     "hatch.linewidth": 0.45, "svg.hashsalt": "HalC8-S4S5"})
INK, MUTED = "#2C2C2A", "#5F5E5A"
FS, FS_S = 7.0, 6.5
W = CFG["width_mm"]
LABEL_W, RIGHT = 43.0, 2.5           # mm
X0, X1 = LABEL_W, W - RIGHT          # map area
KB_L, KB_R = -15.3, 16.3
TRACK_H, ROW_STEP, TOP = 3.0, 15.5, 11.5
GENE_ITALIC = {"proC8", "halU", "halR", "halT1", "halT2", "halT3"}


def X(kb):
    return X0 + (kb - KB_L) / (KB_R - KB_L) * (X1 - X0)


def rich(ax, x, y, segs, size=FS, color=INK, ha="left"):
    """Draw mixed-style text segments [(text, style, weight)] left to right; returns end x."""
    r = ax.figure.canvas.get_renderer()
    inv = ax.transData.inverted()
    widths, arts = [], []
    for t, style, weight in segs:
        a = ax.text(x, y, t, fontsize=size, fontstyle=style, fontweight=weight, color=color, ha="left", va="baseline")
        bb = a.get_window_extent(r)
        widths.append(inv.transform((bb.x1, 0))[0] - inv.transform((bb.x0, 0))[0]); arts.append(a)
    cx = x - (sum(widths) if ha == "right" else 0)
    for a, w in zip(arts, widths):
        a.set_x(cx); cx += w
    return cx


def organism_segs(org):
    parts = org.split(" ")
    segs = [(parts[0], "italic", "normal")]
    rest = parts[1:]
    if rest and rest[0] != "sp.":
        segs.append((" " + rest[0], "italic", "normal")); rest = rest[1:]
    if rest:
        segs.append((" " + " ".join(rest), "normal", "normal"))
    return segs


def arrow(x1, x2, yc, strand, p_left, p_right):
    """Gene arrow polygon between display x1<x2 (mm); jagged edge on a partial end."""
    h = TRACK_H / 2
    head = min(1.4, (x2 - x1) * 0.6)
    def flat(x, top_first, partial):
        ys = [yc + h, yc - h] if top_first else [yc - h, yc + h]
        if not partial:
            return [(x, ys[0]), (x, ys[1])]
        d = min(0.45, (x2 - x1) * 0.3) * (1 if x == x1 else -1)
        n = 4
        return [(x + (d if i % 2 else 0), ys[0] + (ys[1] - ys[0]) * i / n) for i in range(n + 1)]
    if strand > 0:
        left = flat(x1, True, p_left)                     # top -> bottom on the left
        if p_right:
            right = flat(x2, False, True)                 # bottom -> top on the right
            return left + right
        return left + [(x2 - head, yc - h), (x2, yc), (x2 - head, yc + h)]
    right = flat(x2, False, p_right)                      # bottom -> top on the right
    if p_left:
        return right + flat(x1, True, True)
    return right + [(x1 + head, yc + h), (x1, yc), (x1 + head, yc - h)]


def ribbon(ax, a, b, ya, yb, fc, alpha):
    (a1, a2), (b1, b2) = a, b
    ym = (ya + yb) / 2
    verts = [(a1, ya), (a1, ym), (b1, ym), (b1, yb), (b2, yb), (b2, ym), (a2, ym), (a2, ya), (a1, ya)]
    codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MPath(verts, codes), fc=fc, ec="none", alpha=alpha, zorder=1, gid="link"))


def draw(panel):
    rows = defaultdict(list)
    for r in csv.DictReader(open(PKG / "annotations" / f"Figure_{panel}_cds_annotation.csv")):
        for k in ("row", "gbk_start", "gbk_end", "strand", "partial_5", "partial_3"):
            r[k] = int(r[k])
        rows[r["row"]].append(r)
    links = [l for l in csv.DictReader(open(PKG / "annotations" / f"Figure_{panel}_links.tsv"), delimiter="\t") if l["drawn"] == "1"]
    n = len(rows)
    y_axis = TOP + (n - 1) * ROW_STEP + 9.5
    H = y_axis + (31.0 if panel == "S5" else 33.5)
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
    fig.canvas.draw()
    geom, plotted, flipped_any = {}, [], False
    yc_of = {}
    for i in sorted(rows):
        cds = rows[i]
        ref = cds[0]["cluster"] == "REF"
        anchor = next(c for c in cds if (c["cds_id"] == "AAQ82548.1" if ref else c["cds_id"] == c["candidate_accession"]))
        flip = anchor["strand"] < 0
        flipped_any |= flip
        yc = TOP + (i - 1) * ROW_STEP
        yc_of[i] = yc
        def kb(p, _a=anchor, _f=flip):          # p = boundary coordinate (0-based)
            return (_a["gbk_end"] - p) / 1000 if _f else (p - (_a["gbk_start"] - 1)) / 1000
        lo = min(kb(c["gbk_start"] - 1) for c in cds); hi = max(kb(c["gbk_end"]) for c in cds)
        lo2 = min(lo, min(kb(c["gbk_end"]) for c in cds)); hi2 = max(hi, max(kb(c["gbk_start"] - 1) for c in cds))
        ax.add_line(Line2D([X(lo2) - 0.6, X(hi2) + 0.6], [yc, yc], color="#3D3D3A", lw=0.5, zorder=2))
        for c in cds:
            a, b = sorted((kb(c["gbk_start"] - 1), kb(c["gbk_end"])))
            strand = -c["strand"] if flip else c["strand"]
            # partial flags refer to GenBank left (<) / right (>) boundaries
            pl, pr = (c["partial_5"], c["partial_3"]) if not flip else (c["partial_3"], c["partial_5"])
            st = PAL[c["display_class"]]
            ax.add_patch(Polygon(arrow(X(a), X(b), yc, strand, pl, pr), closed=True, fc=st["fill"], ec=st["edge"],
                                 lw=0.45, hatch=st["hatch"] or None, joinstyle="miter", zorder=3,
                                 gid=f"cds|{panel}|{i}|{c['cds_id']}|{c['display_class']}"))
            geom[(i, c["cds_id"])] = (X(a), X(b), c["display_class"])
            plotted.append(dict(c, display_start_kb=f"{a:.3f}", display_end_kb=f"{b:.3f}", display_strand=strand, reverse_complemented=int(flip)))
            if ref:
                up = c["display_class"] in ("halR", "halT1", "halT3")
                xm = (X(a) + X(b)) / 2
                if up:
                    ax.add_line(Line2D([xm, xm], [yc - TRACK_H / 2 - 0.5, yc - TRACK_H / 2 - 3.4], color=MUTED, lw=0.35))
                ax.text(xm, yc - TRACK_H / 2 - (3.9 if up else 1.0), c["display_class"], fontsize=FS, fontstyle="italic",
                        ha="center", va="baseline", color=INK)
        # candidate start guide
        ax.add_line(Line2D([X(0), X(0)], [yc - TRACK_H / 2 - 0.9, yc + TRACK_H / 2 + 0.9], color=INK, lw=0.5, ls=(0, (2, 1.5)), zorder=4))
        # labels
        c0 = cds[0]
        if ref:
            rich(ax, 2.0, yc - 1.2, [("AS7092 reference", "normal", "bold")])
            rich(ax, 2.0, yc + 1.6, organism_segs(c0["organism"]))
            rich(ax, 2.0, yc + 4.3, [("DQ137141.1", "normal", "normal")], size=FS_S, color=MUTED)
        else:
            rich(ax, 2.0, yc - 1.2, [(c0["candidate_accession"] + (" †" if flip else ""), "normal", "bold")])
            rich(ax, 2.0, yc + 1.6, organism_segs(c0["organism"]))
            rich(ax, 2.0, yc + 4.3, [(f"{c0['cluster']} · {c0['nuccore']}", "normal", "normal")], size=FS_S, color=MUTED)
            note = CFG[panel]["row_notes"].get(c0["candidate_accession"])
            if note:
                rich(ax, 2.0, yc + 7.0, [(note, "normal", "normal")], size=FS_S, color="#7A2E14")
    # ribbons
    for l in links:
        a = geom[(int(l["row_a"]), l["cds_a"])]; b = geom[(int(l["row_b"]), l["cds_b"])]
        idn = float(l["identity"])
        fc = PAL[a[2]]["fill"] if a[2] == b[2] and a[2] != "other" else "#8A8985"
        ribbon(ax, a[:2], b[:2], yc_of[int(l["row_a"])] + TRACK_H / 2 + 0.25, yc_of[int(l["row_b"])] - TRACK_H / 2 - 0.25,
               fc, 0.10 + 0.50 * idn)
    # axis
    ax.add_line(Line2D([X(-15), X(16)], [y_axis, y_axis], color=INK, lw=0.6))
    for t in range(-15, 17, 5):
        ax.add_line(Line2D([X(t), X(t)], [y_axis, y_axis + 1.0], color=INK, lw=0.6))
        ax.text(X(t), y_axis + 3.6, ("−" if t < 0 else "+" if t > 0 else "") + str(abs(t)), fontsize=FS, ha="center", va="baseline", color=INK)
    for t in range(-15, 17):
        ax.add_line(Line2D([X(t), X(t)], [y_axis, y_axis + 0.55], color=INK, lw=0.4))
    ax.text((X(-15) + X(16)) / 2, y_axis + 7.0, "Distance from candidate CDS start (kb)", fontsize=FS, ha="center", va="baseline", color=INK)
    # key: three-column grid
    ky, dy, cols = y_axis + 12.5, 5.2, (2.0, 60.0, 116.0)
    def glyph(x, y, cls, partial=0):
        st = PAL[cls]
        ax.add_patch(Polygon(arrow(x, x + 6.0, y, 1, 0, partial), closed=True, fc=st["fill"], ec=st["edge"], lw=0.45,
                             hatch=st["hatch"] or None, zorder=3, gid=f"key|{cls}"))
    for k, cls in enumerate(["halU", "halR", "proC8", "halT1", "halT2", "halT3"]):
        x, y = cols[k % 3], ky + dy * (k // 3)
        glyph(x, y, cls)
        segs = [(cls, "italic", "normal"), (" / ", "normal", "normal"), (cls + "-like", "italic", "normal")]
        if cls == "proC8":
            segs.append((" (candidate)", "normal", "normal"))
        rich(ax, x + 7.2, y + 0.9, segs)
    y = ky + 2 * dy
    if panel == "S5":
        glyph(cols[0], y, "recurrent_neighbour"); rich(ax, cols[0] + 7.2, y + 0.9, [("Recurrent Cluster 4 neighbouring CDS", "normal", "normal")])
        glyph(cols[1], y, "other"); rich(ax, cols[1] + 7.2, y + 0.9, [("Other CDS", "normal", "normal")])
        glyph(cols[2], y, "other", 1); rich(ax, cols[2] + 7.2, y + 0.9, [("Partial CDS", "normal", "normal")])
        y += dy; slots = list(cols)
    else:
        glyph(cols[0], y, "other"); rich(ax, cols[0] + 7.2, y + 0.9, [("Other CDS", "normal", "normal")])
        glyph(cols[1], y, "other", 1); rich(ax, cols[1] + 7.2, y + 0.9, [("Partial CDS", "normal", "normal")])
        slots = [cols[2], cols[0], cols[1]]
    def slot_y(k):
        return y if (panel == "S5" or k == 0) else y + dy
    x = slots[0]
    for j, idn in enumerate((0.3, 0.65, 1.0)):
        ax.add_patch(Rectangle((x + j * 2.0, slot_y(0) - 1.5), 2.0, 3.0, fc="#8A8985", ec="none", alpha=0.10 + 0.50 * idn, gid="key|link"))
    rich(ax, x + 7.2, slot_y(0) + 0.9, [("Clinker link; shade: identity 0.3–1.0", "normal", "normal")])
    x = slots[1]
    ax.add_line(Line2D([x + 3, x + 3], [slot_y(1) - 1.8, slot_y(1) + 1.8], color=INK, lw=0.5, ls=(0, (2, 1.5))))
    rich(ax, x + 7.2, slot_y(1) + 0.9, [("Candidate CDS start (reference: ", "normal", "normal"), ("proC8", "italic", "normal"), (")", "normal", "normal")])
    if flipped_any:
        rich(ax, slots[2], slot_y(2) + 0.9, [("† Locus drawn reverse-complemented", "normal", "normal")])
    out = PKG / "figures"
    meta = {"Creator": "03_plot_figures.py", "CreationDate": None, "ModDate": None}
    fig.savefig(out / f"Figure_{panel}.pdf", metadata=meta)
    fig.savefig(out / f"Figure_{panel}.svg", metadata={"Date": None, "Creator": "03_plot_figures.py"})
    fig.savefig(out / f"Figure_{panel}.png", dpi=600, metadata={"Software": "03_plot_figures.py"})
    plt.close(fig)
    with open(PKG / "annotations" / f"Figure_{panel}_plotted_features.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(plotted[0]), delimiter="\t"); w.writeheader(); w.writerows(plotted)
    print(panel, f"{W}x{H:.1f} mm", len(plotted), "CDS", len(links), "links drawn")


for p in (sys.argv[1:] or ["S4", "S5"]):
    draw(p)
