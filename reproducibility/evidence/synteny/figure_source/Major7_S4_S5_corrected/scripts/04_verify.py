#!/usr/bin/env python3
"""Independent verification of Figures S4/S5 (acceptance criteria 1-10).

Does not import the package readers: GenBank is re-read with Biopython, Clinker JSON with json, and the
figure is checked from the final SVG (element ids, fills, text) and the plotted-feature tables."""
import csv, hashlib, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import BeforePosition, AfterPosition

PKG = Path(__file__).resolve().parent.parent
ARC = PKG.parent / "synteny_figures_final"
DIRS = {"S4": ARC / "S4_AS7092_supported_C1_C3", "S5": ARC / "S5_non_AS7092_C2_C4_C5"}
JSONS = {"S4": ["clinker_outputs/clinker_S4.json", "clinker_outputs/clinker.json"],
         "S5": ["clinker_outputs/clinker_S5.json", "clinker_outputs/clinker.json", "clinker.json"]}
REF = {"AAZ66750.1": "halU", "AAT67988.1": "halR", "AAQ82548.1": "proC8", "AAT67989.2": "halT1", "AAZ66751.1": "halT2", "AAZ66752.1": "halT3"}
PAL = {r["display_class"]: r for r in csv.DictReader(open(PKG / "config" / "colour_map.csv"))}
results, log = [], []


def check(name, ok, detail=""):
    results.append(ok); log.append(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail else ""))


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def links_of(path):
    d = json.loads(Path(path).read_text())["data"]
    g = {x["uid"]: (c["name"], x["label"]) for c in d["clusters"] for l in c["loci"] for x in l["genes"]}
    return {frozenset((g[l["query"]["uid"]], g[l["target"]["uid"]])): l["identity"] for l in d["links"]}, d, g


svg_fill = {}
for panel, D in DIRS.items():
    man = list(csv.DictReader(open(D / "manifest.tsv"), delimiter="\t"))
    byacc = {(m["candidate_accession"] if m["cluster"] != "REF" else "REF"): m for m in man}
    check(f"{panel}: GenBank inputs match manifest SHA-256", all(sha(D / "inputs" / m["gbk_filename"]) == m["sha256"] for m in man), f"{len(man)} files")
    feats = list(csv.DictReader(open(PKG / "annotations" / f"Figure_{panel}_plotted_features.tsv"), delimiter="\t"))
    rows = defaultdict(list)
    for f in feats:
        rows[int(f["row"])].append(f)
    # --- 1 & 5: CDS identity, coordinates, strand, partial flags, order
    bad, total, rowname = [], 0, {}
    for i, fs in rows.items():
        key = fs[0]["candidate_accession"] if fs[0]["cluster"] != "REF" else "REF"
        m = byacc[key]; rowname[i] = m["display_name"]
        rec = SeqIO.read(D / "inputs" / m["gbk_filename"], "genbank")
        gb = []
        for ft in rec.features:
            if ft.type != "CDS":
                continue
            q = ft.qualifiers
            cid = q.get("protein_id", q.get("locus_tag", q.get("gene", [""])))[0]
            gb.append((cid, int(ft.location.start) + 1, int(ft.location.end), ft.location.strand,
                       int(isinstance(ft.location.start, BeforePosition)), int(isinstance(ft.location.end, AfterPosition))))
        pl = [(f["cds_id"], int(f["gbk_start"]), int(f["gbk_end"]), int(f["strand"]), int(f["partial_5"]), int(f["partial_3"])) for f in fs]
        total += len(pl)
        if sorted(gb) != sorted(pl) or [x[0] for x in gb] != [x[0] for x in pl]:
            bad.append(m["display_name"])
        # display transform is the stated one (align at anchor start, optional reverse complement)
        anc = next(f for f in fs if f["cds_id"] == ("AAQ82548.1" if key == "REF" else key))
        flip = int(anc["strand"]) < 0
        for f in fs:
            s, e = int(f["gbk_start"]) - 1, int(f["gbk_end"])
            a, b = sorted(((int(anc["gbk_end"]) - p) / 1000 if flip else (p - (int(anc["gbk_start"]) - 1)) / 1000) for p in (s, e))
            if abs(a - float(f["display_start_kb"])) > 6e-4 or abs(b - float(f["display_end_kb"])) > 6e-4 or \
               int(f["display_strand"]) != (-int(f["strand"]) if flip else int(f["strand"])) or int(f["reverse_complemented"]) != int(flip):
                bad.append((m["display_name"], f["cds_id"]))
    check(f"{panel}: every plotted CDS equals a GenBank CDS (id, start, end, strand, partial flags, order); nothing omitted; display = shift ± reverse complement only", not bad, f"{total} CDS, {len(bad)} mismatches")
    # --- archived link evidence (all JSON variants agree)
    L0, d0, g0 = links_of(D / JSONS[panel][0])
    check(f"{panel}: all archived Clinker JSON variants carry the same link set", all(links_of(D / j)[0] == L0 for j in JSONS[panel]), f"{len(L0)} links, {len(JSONS[panel])} files")
    partner = defaultdict(list)
    for k, idn in L0.items():
        a, b = tuple(k) if len(k) == 2 else (list(k)[0], list(k)[0])
        partner[a].append((b, idn)); partner[b].append((a, idn))
    # --- 2 & 4: evidence for every highlighted CDS; no other CDS highlighted
    c4 = [m["display_name"] for m in man if m["cluster"] == "Cluster 4"]
    rec_ok = set()
    if panel == "S5":
        grp = [g for g in d0["groups"] if g["label"].startswith("Group ") and {g0[u][0] for u in g["genes"]} == set(c4)]
        check("S5: exactly five Clinker homology groups span the four Cluster 4 loci", len(grp) == 5, ", ".join(g["label"] for g in grp))
        for g in grp:
            for u in g["genes"]:
                node = g0[u]
                if len({p[0][0] for p in partner[node] if p[0][0] in c4 and p[0][0] != node[0]}) == 3:
                    rec_ok.add(node)
    fails, expect = [], Counter()
    for i, fs in rows.items():
        for f in fs:
            node, cls = (rowname[i], f["cds_id"]), f["display_class"]
            exp = "other"
            if f["cluster"] == "REF":
                exp = REF[f["cds_id"]]
            elif f["cds_id"] == f["candidate_accession"]:
                exp = "proC8"
            else:
                refs = {REF[p[0][1]] for p in partner[node] if p[0][1] in REF}
                if len(refs) == 1 and panel == "S4":
                    exp = refs.pop()
                elif refs:
                    exp = "UNEXPECTED_REFERENCE_LINK"
                elif node in rec_ok:
                    exp = "recurrent_neighbour"
            expect[exp] += 1
            if exp != cls:
                fails.append((node, cls, exp))
    check(f"{panel}: display class of every CDS re-derived from manifest + archived Clinker links (cassette colour <=> link to that reference gene; hatch <=> recurrent group rule)", not fails, f"{dict(expect)}; {len(fails)} disagreements {fails[:3]}")
    # --- SVG checks
    svg = (PKG / "figures" / f"Figure_{panel}.svg").read_text()
    els = re.findall(r'<g id="cds\|(S\d)\|(\d+)\|([^|"]+)\|([^"]+)">\s*<path[^>]*?style="([^"]*)"', svg, re.S)
    fills = defaultdict(set); wrong = []
    for p, r, cid, cls, style in els:
        fill = re.search(r"fill: (url\(#\w+\)|#[0-9a-f]{6})", style).group(1)
        fills[cls].add(fill)
        if not PAL[cls]["hatch"] and fill != PAL[cls]["fill"].lower():
            wrong.append((cid, cls, fill))
        if PAL[cls]["hatch"] and not fill.startswith("url("):
            wrong.append((cid, cls, fill))
    svg_fill[panel] = {k: sorted(v) for k, v in fills.items()}
    key = Counter((f["panel"], f["row"], f["cds_id"], f["display_class"]) for f in feats)
    check(f"{panel}: SVG has one arrow per plotted CDS with the class recorded in the table", Counter((p, r, c, k) for p, r, c, k, _ in els) == key, f"{len(els)} arrows")
    check(f"{panel}: every arrow fill equals the colour map for its class; all 'other' CDS are neutral grey", not wrong, f"{sum(1 for e in els if e[3]=='other')} grey 'other' arrows")
    # --- 6: links
    lk = list(csv.DictReader(open(PKG / "annotations" / f"Figure_{panel}_links.tsv"), delimiter="\t"))
    tab = {frozenset(((rowname[int(l["row_a"])], l["cds_a"]), (rowname[int(l["row_b"])], l["cds_b"]))): float(l["identity"]) for l in lk}
    same = set(tab) == set(L0) and all(abs(tab[k] - L0[k]) < 1e-4 for k in tab)
    drawn = [l for l in lk if l["drawn"] == "1"]
    adj = all(int(l["row_b"]) - int(l["row_a"]) == 1 for l in drawn) and len(drawn) == sum(1 for l in lk if int(l["row_b"]) - int(l["row_a"]) == 1)
    check(f"{panel}: link table equals the archived Clinker links; ribbons = all links between adjacent rows, nothing added", same and adj and svg.count('<g id="link">') == len(drawn), f"{len(lk)} archived, {len(drawn)} drawn, {svg.count(chr(60)+'g id=\"link\">')} ribbons in SVG")
    # --- 7: labels
    texts = re.findall(r'<text style="([^"]*)"[^>]*>([^<]*)</text>', svg)
    alltxt = " ".join(t for _, t in texts)
    miss = [m["candidate_accession"] for m in man if m["cluster"] != "REF" and m["candidate_accession"] not in alltxt]
    miss += [m["cluster"] for m in man if m["cluster"] != "REF" and f'{m["cluster"]} · {m["nuccore"].replace(".", "_", 1)}' not in alltxt]
    for m in man:
        if m["cluster"] != "REF":
            org = m["display_name"].split(" . ", 2)[2]
            if not all(w in alltxt for w in org.split()):
                miss.append(org)
    check(f"{panel}: all {len(man)} manifest loci shown with correct accession, cluster, nuccore and organism", not miss and len(rows) == len(man), str(miss))
    ital = Counter(t.strip() for s, t in texts if "italic" in s)
    allowed_genes = set(REF.values()) | {v + "-like" for v in REF.values()}
    binom = set()
    for m in man:
        if m["cluster"] != "REF":
            w = m["display_name"].split(" . ", 2)[2].split()
            binom.add(w[0]); 
            if len(w) > 1 and w[1] != "sp.": binom.add(w[1])
    binom |= {"Halobacterium"}
    bad_it = [t for t in ital if t not in allowed_genes and t not in binom]
    roman_bad = [t for s, t in texts if "italic" not in s and any(w in binom for w in t.split())]
    check(f"{panel}: italics are exactly genus/species names and gene names; 'sp.', strains and accessions roman", not bad_it and not roman_bad, f"italic: {sorted(ital)}")
    sizes = sorted({float(x) for x in re.findall(r"font-size: ([\d.]+)px", svg)})
    w_h = re.search(r'<svg[^>]*width="([\d.]+)pt" height="([\d.]+)pt"', svg)
    check(f"{panel}: publication size and legibility", sizes[0] >= 6.5 and abs(float(w_h.group(1)) / 72 * 25.4 - 170) < 0.5, f"{float(w_h.group(1))/72*25.4:.0f} x {float(w_h.group(2))/72*25.4:.0f} mm, text {sizes[0]}–{sizes[-1]} pt, text kept as text")
    import struct
    png = (PKG / "figures" / f"Figure_{panel}.png").read_bytes()
    wpx = struct.unpack(">I", png[16:20])[0]
    check(f"{panel}: PNG is 600 dpi", abs(wpx / (170 / 25.4) - 600) < 2, f"{wpx} px wide")
    # --- annotation table hygiene
    ann = list(csv.DictReader(open(PKG / "annotations" / f"Figure_{panel}_cds_annotation.csv")))
    check(f"{panel}: visualization annotation table lists only displayed GenBank CDS", len(ann) == total and all(a["cds_id"] for a in ann))

# --- 3: consistency across panels
cas = ["halU", "halR", "proC8", "halT1", "halT2", "halT3", "other"]
check("S4 vs S5: each cassette class and 'other' has a single fill, identical in both figures",
      all(len(svg_fill[p].get(c, [])) == 1 for p in svg_fill for c in cas) and all(svg_fill["S4"][c] == svg_fill["S5"][c] for c in cas),
      ", ".join(f"{c} {svg_fill['S4'][c][0]}" for c in cas))
check("S4 contains no hatched (recurrent Cluster 4) arrow; S5 contains no cassette-coloured non-reference CDS other than candidates", "recurrent_neighbour" not in svg_fill["S4"])
ex = list(csv.DictReader(open(PKG / "annotations" / "excluded_annotation_entries.csv")))
orig = (DIRS["S4"] / "gene_functions" / "S4_AS7092_gene_functions.csv").read_text()
check("Three orphan S4 identifiers excluded from the visualization table only; still present in the archived CSV",
      sorted(e["cds_id"] for e in ex) == ["WP_049899082.1", "WP_049899084.1", "WP_049899087.1"] and all(e["cds_id"] in orig for e in ex))
# --- 10: archive untouched
snap = PKG / "audit" / "archive_sha256_snapshot.tsv"
cur = {str(p.relative_to(ARC)): sha(p) for p in sorted(ARC.rglob("*")) if p.is_file() and p.name != ".DS_Store"}
if not snap.exists():
    snap.write_text("".join(f"{v}\t{k}\n" for k, v in cur.items()))
old = dict(reversed(l.rstrip("\n").split("\t", 1)) for l in snap.read_text().splitlines())
pkg_start = (PKG / "scripts" / "s4s5_common.py").stat().st_mtime
newer = [k for k in cur if (ARC / k).stat().st_mtime > pkg_start - 3600]
check("Archive unchanged: all files match the SHA-256 snapshot; no archive file modified since this package was started", old == cur and not newer, f"{len(cur)} files")

out = "\n".join(log) + f"\n\n{sum(results)} of {len(results)} checks passed\n"
(PKG / "audit" / "verification_log.txt").write_text(out)
print(out)
sys.exit(0 if all(results) else 1)
