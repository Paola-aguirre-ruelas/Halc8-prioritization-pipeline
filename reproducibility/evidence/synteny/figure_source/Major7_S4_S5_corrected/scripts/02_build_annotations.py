#!/usr/bin/env python3
"""Build the visualization-specific annotation tables for Figures S4 and S5.

Every displayed CDS comes from an archived GenBank window. A CDS is highlighted only when the
archive documents why (rules below); everything else is 'other'. Archive files are read only.

  reference gene      /gene qualifier of the AS7092 reference record (DQ137141.1)
  proC8 (candidate)   protein_id equals the manifest candidate_accession
  <gene>-like         labelled so in the archived panel annotation CSV, present in the displayed
                      GenBank window, and joined to that reference gene by an archived Clinker link
  recurrent_neighbour S5 only: member of a Clinker homology group found in all four Cluster 4 loci,
                      and itself linked to a member in each of the three other Cluster 4 loci
"""
import csv, json
from collections import defaultdict
from s4s5_common import *

OUT = PKG / "annotations"
CFG = json.loads((PKG / "config" / "figure_config.json").read_text())
FIELDS = ["panel", "row", "cluster", "candidate_accession", "nuccore", "organism", "cds_id", "protein_id",
          "locus_tag", "gbk_start", "gbk_end", "strand", "partial_5", "partial_3", "pseudo", "product",
          "archived_label", "display_class", "evidence"]
excluded = []

for panel, P in PANELS.items():
    man = {(r["candidate_accession"] if r["cluster"] != "REF" else "REF"): r for r in read_manifest(panel)}
    fx = dict(read_functions(panel))
    raw = load_clinker(P["dir"] / P["raw_json"])
    G = raw["genes"]
    # links indexed by (cluster name, label)
    nb = defaultdict(list)
    for l in raw["links"]:
        q, t = G[l["q"]], G[l["t"]]
        nb[(q["cluster"], q["label"])].append((t["cluster"], t["label"], l["identity"], l["similarity"]))
        nb[(t["cluster"], t["label"])].append((q["cluster"], q["label"], l["identity"], l["similarity"]))
    c4 = [r["display_name"] for r in man.values() if r["cluster"] == "Cluster 4"]
    rec_groups = {}
    if len(c4) == 4:
        for g in raw["groups"]:
            if not g["label"].startswith("Group "):
                continue
            if {G[u]["cluster"] for u in g["genes"]} == set(c4):
                for u in g["genes"]:
                    rec_groups[(G[u]["cluster"], G[u]["label"])] = g["label"]
        assert len(set(rec_groups.values())) == 5, sorted(set(rec_groups.values()))
    rows, shown_ids = [], set()
    for i, key in enumerate(CFG[panel]["row_order"], 1):
        r = man[key]
        rec = parse_gbk(P["dir"] / "inputs" / r["gbk_filename"])
        name = r["display_name"]
        org = "Halobacterium sp. AS7092" if key == "REF" else name.split(" . ", 2)[2]
        for c in rec["cds"]:
            cid = cds_id(c)
            shown_ids.add(cid)
            links = nb.get((name, cid), [])
            arch = fx.get(cid, "")
            cls, ev = "other", "no highlight criterion met"
            if key == "REF":
                cls = REF_GENES[cid]
                ev = f"GenBank /gene={c['gene']} in the AS7092 reference record (DQ137141.1)"
                assert c["gene"] == cls and arch == cls
            elif cid == r["candidate_accession"]:
                cls = "proC8"
                ref = [x for x in links if x[1] == "AAQ82548.1"]
                ev = "manifest candidate_accession; " + (
                    f"Clinker link to proC8 AAQ82548.1, identity {ref[0][2]:.2f}" if ref else
                    "no archived Clinker link to proC8 AAQ82548.1 (candidate status comes from the scoring pipeline)")
                assert arch == "proC8-like"
            elif arch.endswith("-like"):
                want = arch[:-5]
                ref = [x for x in links if REF_GENES.get(x[1]) == want]
                assert ref, (cid, arch)
                cls = want
                ev = f"archived label {arch}; Clinker link to {want} {ref[0][1]}, identity {ref[0][2]:.2f}"
            elif (name, cid) in rec_groups:
                others = {x[0] for x in links if x[0] in c4 and x[0] != name}
                grp = rec_groups[(name, cid)]
                if len(others) == 3:
                    cls = "recurrent_neighbour"
                    mn = min(x[2] for x in links if x[0] in c4)
                    ev = f"Clinker {grp}; linked to a CDS in each of the 3 other Cluster 4 loci (minimum identity {mn:.2f})"
                else:
                    ev = (f"member of Clinker {grp} but linked to only {len(others)} of 3 other Cluster 4 loci "
                          f"(max identity {max(x[2] for x in links):.2f}); left unhighlighted")
            rows.append(dict(panel=panel, row=i, cluster=r["cluster"], candidate_accession=r["candidate_accession"],
                             nuccore=r["nuccore"].replace(".", "_", 1) if key != "REF" else r["nuccore"],
                             organism=org, cds_id=cid, protein_id=c["protein_id"], locus_tag=c["locus_tag"],
                             gbk_start=c["start"], gbk_end=c["end"], strand=c["strand"],
                             partial_5=int(c["partial5"]), partial_3=int(c["partial3"]), pseudo=int(c["pseudo"]),
                             product=c["product"], archived_label=arch or "(not listed)", display_class=cls, evidence=ev))
    with open(OUT / f"Figure_{panel}_cds_annotation.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
    for cid, lab in read_functions(panel):
        if cid not in shown_ids:
            excluded.append(dict(panel=panel, cds_id=cid, archived_label=lab,
                                 reason="identifier absent from every displayed GenBank window of this panel; "
                                        "kept in the archived CSV, omitted from the visualization table only"))
    # all archived links between displayed loci, with drawn/not-drawn flag (adjacent rows only)
    order = [man[k]["display_name"] for k in CFG[panel]["row_order"]]
    with open(OUT / f"Figure_{panel}_links.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["panel", "row_a", "cds_a", "row_b", "cds_b", "identity", "similarity", "drawn"])
        for l in raw["links"]:
            q, t = G[l["q"]], G[l["t"]]
            a, b = order.index(q["cluster"]) + 1, order.index(t["cluster"]) + 1
            if a > b:
                a, b, q, t = b, a, t, q
            w.writerow([panel, a, q["label"], b, t["label"], f"{l['identity']:.4f}", f"{l['similarity']:.4f}",
                        int(b - a == 1 and l["identity"] >= CFG["min_link_identity"])])
    from collections import Counter
    print(panel, len(rows), "CDS;", dict(Counter(r["display_class"] for r in rows)))
    for r in rows:
        if "left unhighlighted" in r["evidence"]:
            print("   note:", r["cds_id"], r["evidence"])

with open(OUT / "excluded_annotation_entries.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["panel", "cds_id", "archived_label", "reason"]); w.writeheader(); w.writerows(excluded)
print("excluded:", [(e["panel"], e["cds_id"]) for e in excluded])
import shutil
shutil.copy(PKG / "config" / "colour_map.csv", OUT / "colour_map_S4_S5.csv")
