#!/usr/bin/env python3
"""Step 1-2: source integrity and independent reproduction of the reported graphical defects.
Reads the archive only; writes audit/source_audit.json and audit/old_svg_gene_colours_<panel>.tsv."""
import csv, json, re, sys
from collections import Counter, defaultdict
from s4s5_common import *

OUT = PKG / "audit"
GREY = {"#b0b0b0", "#bbb", "#bbbbbb"}


def norm_col(c):
    c = c.strip().lower()
    m = re.match(r"rgb\((\d+),\s*(\d+),\s*(\d+)\)", c)
    if m:
        return "#%02x%02x%02x" % tuple(int(x) for x in m.groups())
    return c


def is_red(hexc):
    if not re.fullmatch(r"#[0-9a-f]{6}", hexc):
        return False
    r, g, b = (int(hexc[i:i + 2], 16) for i in (1, 3, 5))
    return r >= 180 and g <= 60 and b <= 60


def svg_genes(path):
    """Gene polygons in document order: (clinker uid, fill, group uid, label text, row label)."""
    s = open(path).read()
    out, row = [], ""
    tok = re.compile(r'<text class="clusterText"[^>]*>([^<]*)<|<g id="gene_([0-9a-f-]+)" class="gene"[^>]*>\s*'
                     r'<polygon class="genePolygon([^"]*)"([^>]*)>\s*<text class="geneLabel"[^>]*>([^<]*)<')
    for m in tok.finditer(s):
        if m.group(1) is not None:
            row = m.group(1).strip()
            continue
        attrs = m.group(4)
        fill = re.search(r'fill="([^"]+)"', attrs) or re.search(r'fill:\s*([^;"]+)', attrs)
        grp = re.search(r"group-([0-9a-f-]+)", m.group(3))
        out.append((m.group(2), norm_col(fill.group(1)) if fill else "none", grp.group(1) if grp else "", m.group(5).strip(), row))
    return out


report = {}
for panel, P in PANELS.items():
    R = report[panel] = {}
    man = read_manifest(panel)
    # 1. checksums
    R["checksums"] = [dict(file=r["gbk_filename"], ok=sha256(P["dir"] / "inputs" / r["gbk_filename"]) == r["sha256"]) for r in man]
    # GenBank CDS
    gbk = {r["display_name"]: parse_gbk(P["dir"] / "inputs" / r["gbk_filename"]) for r in man}
    gb_ids = defaultdict(list)
    for name, rec in gbk.items():
        for c in rec["cds"]:
            gb_ids[cds_id(c)].append(name)
    R["gbk_cds_per_locus"] = {n: len(r["cds"]) for n, r in gbk.items()}
    R["gbk_cds_without_protein_id"] = {n: [c["locus_tag"] for c in r["cds"] if not c["protein_id"]] for n, r in gbk.items()}
    # Clinker JSON variants
    raw = load_clinker(P["dir"] / P["raw_json"])
    def linkset(j):
        return sorted((j["genes"][l["q"]]["label"], j["genes"][l["t"]]["label"], round(l["identity"], 6)) for l in j["links"])
    R["json_variants"] = []
    for o in [P["raw_json"]] + P["other_json"]:
        j = load_clinker(P["dir"] / o)
        R["json_variants"].append(dict(
            file=o, sha256=sha256(P["dir"] / o), n_genes=len(j["genes"]), n_links=len(j["links"]),
            same_links_as_raw=linkset(j) == linkset(raw),
            same_gene_coords_after_unflip=all(
                (lambda g, L: ((L["end"] - g["end"], L["end"] - g["start"], -g["strand"]) if L["flipped"] else (g["start"], g["end"], g["strand"])))(g, j["loci"][g["locus"]])
                == (raw["genes"][u]["start"], raw["genes"][u]["end"], raw["genes"][u]["strand"]) for u, g in j["genes"].items()),
            flipped_loci=[v["cluster"] for v in j["loci"].values() if v["flipped"]],
            slot_order=[n for n, s in sorted(j["clusters"], key=lambda x: x[1])],
            functional_group_colours={g["label"]: g["colour"] for g in j["groups"] if not g["label"].startswith("Group ")}))
    # GenBank <-> Clinker JSON correspondence (raw JSON: 0-based start, end as GenBank end)
    mism = []
    jl = defaultdict(list)
    for u, g in raw["genes"].items():
        jl[g["cluster"]].append(g)
    for name, rec in gbk.items():
        a = sorted((cds_id(c), c["start"] - 1, c["end"], c["strand"]) for c in rec["cds"])
        b = sorted((g["label"], g["start"], g["end"], g["strand"]) for g in jl[name])
        if a != b:
            mism.append(dict(locus=name, only_gbk=[x for x in a if x not in b][:10], only_json=[x for x in b if x not in a][:10]))
    R["gbk_vs_clinker_json_mismatches"] = mism
    R["clinker_genes_total"] = len(raw["genes"])
    # annotation CSV vs GenBank
    fx = read_functions(panel)
    fmap = dict(fx)
    R["csv_rows"] = len(fx)
    R["csv_duplicates"] = [k for k, v in Counter(i for i, _ in fx).items() if v > 1]
    R["csv_ids_not_in_gbk"] = [(i, f) for i, f in fx if i not in gb_ids]
    R["gbk_ids_not_in_csv"] = sorted(i for i in gb_ids if i not in fmap)
    R["csv_class_counts"] = Counter(f for _, f in fx)
    # cassette labels vs Clinker links to the reference
    best = defaultdict(dict)
    for l in raw["links"]:
        q, t = raw["genes"][l["q"]], raw["genes"][l["t"]]
        for a, b in ((q, t), (t, q)):
            if a["label"] in REF_GENES:
                best[b["label"]][REF_GENES[a["label"]]] = round(l["identity"], 3)
    lab = []
    for i, f in fx:
        if f.endswith("-like") and i in gb_ids:
            want = f[:-5]
            lab.append(dict(id=i, label=f, locus=gb_ids[i][0], links_to_reference=best.get(i, {}), supported=want in best.get(i, {})))
    R["like_labels"] = lab
    R["like_labels_without_reference_link"] = [x for x in lab if not x["supported"]]
    R["other_or_unlisted_CDS_with_reference_link"] = {i: v for i, v in best.items() if fmap.get(i, "other CDS") == "other CDS"}
    # old SVG colours
    R["old_svg"] = {}
    for svg in P["old_svgs"]:
        path = svg if hasattr(svg, "exists") and svg.is_absolute() else P["dir"] / svg
        rows = []
        for uid, fill, grp, label, row in svg_genes(path):
            cls = fmap.get(label, "NOT_IN_CSV")
            rows.append(dict(gene_uid=uid, uid_in_raw_json=uid in raw["genes"], cds=label, locus=row,
                             in_gbk=label in gb_ids, csv_class=cls, fill=fill))
        with open(OUT / f"old_svg_gene_colours_{panel}_{path.name}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
        oth = [r for r in rows if r["csv_class"] in ("other CDS", "NOT_IN_CSV")]
        col = [r for r in oth if r["fill"] not in GREY]
        cas = defaultdict(Counter)
        for r in rows:
            if r["csv_class"] not in ("other CDS", "NOT_IN_CSV"):
                cas[r["csv_class"].replace("-like", "")][r["fill"]] += 1
        R["old_svg"][str(path.relative_to(ARCHIVE))] = dict(
            polygons=len(rows), uids_in_raw_json=sum(r["uid_in_raw_json"] for r in rows), labels_in_gbk=sum(r["in_gbk"] for r in rows),
            other_cds_polygons=len(oth), other_cds_coloured=len(col),
            other_cds_coloured_listed_in_csv=sum(r["csv_class"] == "other CDS" for r in col),
            other_cds_coloured_not_in_csv=sum(r["csv_class"] == "NOT_IN_CSV" for r in col),
            other_cds_red=[(r["cds"], r["locus"], r["fill"]) for r in col if is_red(r["fill"])],
            other_cds_colour_counts=Counter(r["fill"] for r in col).most_common(),
            cassette_fill_by_gene={k: dict(v) for k, v in cas.items()})
(OUT / "source_audit.json").write_text(json.dumps(report, indent=1, default=str))

for panel, R in report.items():
    print("=====", panel)
    print("checksums ok:", sum(x["ok"] for x in R["checksums"]), "/", len(R["checksums"]))
    print("CDS per locus:", R["gbk_cds_per_locus"]); print("no protein_id:", R["gbk_cds_without_protein_id"])
    for v in R["json_variants"]:
        print(" json", v["file"], v["sha256"][:10], v["n_genes"], v["n_links"], "links==raw", v["same_links_as_raw"], "coords", v["same_gene_coords_after_unflip"], "flipped", len(v["flipped_loci"]), v["functional_group_colours"])
    print("gbk vs json mismatches:", R["gbk_vs_clinker_json_mismatches"])
    print("csv rows", R["csv_rows"], "dups", R["csv_duplicates"], "classes", dict(R["csv_class_counts"]))
    print("csv ids not in gbk:", R["csv_ids_not_in_gbk"]); print("gbk ids not in csv:", R["gbk_ids_not_in_csv"])
    print("-like without ref link:", R["like_labels_without_reference_link"])
    print("other CDS with ref link:", R["other_or_unlisted_CDS_with_reference_link"])
    for k, v in R["old_svg"].items():
        print(" svg", k, {a: b for a, b in v.items() if a != "other_cds_colour_counts"}); print("   colours", v["other_cds_colour_counts"][:8])
