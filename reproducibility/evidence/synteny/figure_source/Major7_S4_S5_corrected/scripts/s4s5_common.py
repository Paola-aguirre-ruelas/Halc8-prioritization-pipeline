"""Shared readers for the Major 7 S4/S5 visualization package (read-only on the archive)."""
import csv, hashlib, json, re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PKG = HERE.parent
ARCHIVE = PKG.parent / "synteny_figures_final"
PANELS = {
    "S4": dict(dir=ARCHIVE / "S4_AS7092_supported_C1_C3",
               csv="S4_AS7092_gene_functions.csv",
               raw_json="clinker_outputs/clinker_S4.json",
               other_json=["clinker_outputs/clinker.json"],
               old_svgs=[ARCHIVE / "Fig_S4_PREVIEW.svg", "clinker_outputs/clinker.svg"]),
    "S5": dict(dir=ARCHIVE / "S5_non_AS7092_C2_C4_C5",
               csv="S5_non_AS7092_gene_functions.csv",
               raw_json="clinker_outputs/clinker_S5.json",
               other_json=["clinker_outputs/clinker.json", "clinker.json"],
               old_svgs=[ARCHIVE / "Fig_S5.svg", "clinker.svg"]),
}
REF_GENES = {"AAZ66750.1": "halU", "AAT67988.1": "halR", "AAQ82548.1": "proC8",
             "AAT67989.2": "halT1", "AAZ66751.1": "halT2", "AAZ66752.1": "halT3"}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_manifest(panel):
    with open(PANELS[panel]["dir"] / "manifest.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return sorted(rows, key=lambda r: int(r["sort_order"]))


def read_functions(panel):
    p = PANELS[panel]["dir"] / "gene_functions" / PANELS[panel]["csv"]
    with open(p, newline="") as fh:
        return [(r[0].strip(), r[1].strip()) for r in csv.reader(fh) if r]


_LOC = re.compile(r"(\d+)")


def parse_gbk(path):
    """Minimal GenBank reader: record length and CDS features (1-based inclusive)."""
    text = Path(path).read_text()
    m = re.search(r"^LOCUS\s+(\S+)\s+(\d+) bp", text, re.M)
    rec = dict(locus=m.group(1), length=int(m.group(2)), cds=[])
    m = re.search(r"^ACCESSION\s+(\S+)", text, re.M)
    rec["accession"] = m.group(1) if m else ""
    m = re.search(r"^VERSION\s+(\S+)", text, re.M)
    rec["version"] = m.group(1) if m else ""
    m = re.search(r"^  ORGANISM\s+(.+)$", text, re.M)
    rec["organism"] = m.group(1).strip() if m else ""
    feats = text.split("\nFEATURES", 1)[1].split("\nORIGIN", 1)[0]
    blocks = re.split(r"\n(?= {5}\S)", feats)
    for b in blocks:
        lines = b.split("\n")
        head = lines[0]
        if not head.startswith("     CDS "):
            continue
        loc = head[21:].strip()
        i = 1
        while i < len(lines) and not lines[i].strip().startswith("/"):
            loc += lines[i].strip()
            i += 1
        quals = {}
        body = "\n".join(lines[i:])
        for qm in re.finditer(r'/(\w+)(?:=(?:"((?:[^"]|"")*)"|(\S+)))?', body):
            k = qm.group(1)
            v = qm.group(2) if qm.group(2) is not None else (qm.group(3) or "")
            quals.setdefault(k, re.sub(r"\s*\n\s*", " ", v))
        nums = [int(x) for x in _LOC.findall(loc)]
        rec["cds"].append(dict(
            location=loc, start=min(nums), end=max(nums),
            strand=-1 if "complement" in loc else 1,
            partial5="<" in loc, partial3=">" in loc,
            protein_id=quals.get("protein_id", ""), locus_tag=quals.get("locus_tag", ""),
            gene=quals.get("gene", ""), product=quals.get("product", ""),
            pseudo=("pseudo" in quals or "pseudogene" in quals)))
    return rec


def cds_id(c):
    return c["protein_id"] or c["locus_tag"] or c["gene"]


def load_clinker(path):
    d = json.loads(Path(path).read_text())["data"]
    genes, loci = {}, {}
    for c in d["clusters"]:
        for l in c["loci"]:
            loci[l["uid"]] = dict(cluster=c["name"], flipped=bool(l.get("_flipped")),
                                  start=l["start"], end=l["end"])
            for g in l["genes"]:
                genes[g["uid"]] = dict(label=g["label"], cluster=c["name"], start=g["start"],
                                       end=g["end"], strand=g["strand"], locus=l["uid"])
    links = [dict(uid=l["uid"], q=l["query"]["uid"], t=l["target"]["uid"],
                  identity=l["identity"], similarity=l["similarity"]) for l in d["links"]]
    groups = [dict(uid=g["uid"], label=g["label"], colour=g.get("colour"), genes=g["genes"])
              for g in d["groups"]]
    return dict(genes=genes, loci=loci, links=links, groups=groups,
                clusters=[(c["name"], c.get("slot")) for c in d["clusters"]])
