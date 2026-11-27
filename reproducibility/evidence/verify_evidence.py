#!/usr/bin/env python3
"""Verify archived evidence and rerun selected computations. Does not rerun external bioinformatics tools."""
from pathlib import Path
import csv, hashlib, tempfile, subprocess, sys, importlib.util, collections
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent

def check(ok, label):
    if not ok: raise RuntimeError(label)
    print('PASS: '+label)

def rows(p, delimiter=','):
    with p.open(newline='') as f: return list(csv.DictReader(f, delimiter=delimiter))

for line in (HERE/'SHA256SUMS.txt').read_text().splitlines():
    h,p=line.split('  ',1)
    if hashlib.sha256((HERE/p).read_bytes()).hexdigest()!=h: raise RuntimeError('Hash mismatch: '+p)
print('Evidence integrity: PASS')
e=HERE/'extraction'; spec=importlib.util.spec_from_file_location('extractor', e/'extract_refseq_domains_posthoc.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
with (e/'refseq44_source_proteins.faa').open() as f,(HERE/'retrieval/results/halobacteria_HalC8_v2.2.iEvalue_1e-10.domtblout').open() as d: observed=m.extract(f,d)
with (e/'halobacteria_44_HalC8_domains_extracted_by_hmmer.faa').open() as f: expected=dict(m.fasta_records(f))
with (e/'tabla5_tfm2021_HalC8_domains_extracted_by_hmmer.faa').open() as f: calibration=dict(m.fasta_records(f))
check(observed==expected and len(observed)==44,'44 reconstructed domains')
novel={k:v for k,v in observed.items() if v not in set(calibration.values())}
check(len(novel)==18 and len(calibration)==42,'26 exact calibration matches; 18 retained')
check(novel=={r['accession']:r['sequence'] for r in rows(ROOT/'reproducibility/data/refseq18_for_scoring.csv')},'18 scoring inputs match reconstructed domains')
controls=HERE/'controls/audit_S7_canonical_20260614'
h=controls/'rerun_hard_negative_legacy'; g=controls/'rerun_genome_pilot'; r=HERE/'comparators/audit_S9_de_novo_20260615/01_ri2023'
cases=[(h/'negative_set_legacy_domain_compatible_CANONICAL_input.csv',h/'negative_set_domain_compatible_CANONICAL_scores.csv',8,0),(h/'negative_set_legacy9_all_CANONICAL_input.csv',h/'negative_set_legacy9_all_CANONICAL_scores.csv',9,0),(g/'genome_pilot_length40_120_CANONICAL_input.csv',g/'genome_pilot_length40_120_CANONICAL_scores.csv',2062,4),(r/'ri2023_OP548098_translated_CDS_CANONICAL_input.csv',r/'ri2023_OP548098_translated_CDS_CANONICAL_scores.csv',1,1)]
with tempfile.TemporaryDirectory() as d:
    for i,(inp,expected,n,pos) in enumerate(cases):
        out=Path(d)/f'{i}.csv'
        subprocess.run([sys.executable,str(ROOT/'halc8_score_domains.py'),'--input',str(inp),'--output',str(out),'--input-type','domain','--ss-mode','none'],check=True,stdout=subprocess.DEVNULL)
        result=rows(out)
        check(out.read_bytes()==expected.read_bytes() and len(result)==n and sum(x['verdict_v2'] in ('PROBABLE','MOST LIKELY') for x in result)==pos,expected.name)
cl=HERE/'clustering/clustering_n13_canonical_20260613/tables/cluster_counts_n13_primary.csv'
check([(r['cdhit_clusters'],r['mmseqs_clusters']) for r in rows(cl)]==[(str(i),str(i)) for i in (9,7,5,4)],'archived n13 cluster counts')
loc=rows(HERE/'synteny/reconciliation/locus_inventory_56_RECONCILED.tsv','\t')
check(len(loc)==56 and len({(r['candidate_accession'],r['nuccore_accession']) for r in loc})==56 and sum(r['evaluable_status']=='EVALUABLE_FULL_PM15KB' for r in loc)==42,'56 locus records; 42 full')
# TG-1: verify the published five-record comparison against existing inputs.
t=HERE/'comparators/tg1'
qrows=rows(t/'query.csv'); q=qrows[0]['sequence']
check(len(qrows)==1 and len(q)==77,'TG-1 query length')
map_rows=rows(HERE/'clustering/cluster_id_mapping.tsv','\t')
historical={r['historical_id'] for r in map_rows if r['manuscript_id']=='4'}
members=rows(HERE/'clustering/clustering_n13_canonical_20260613/tables/cluster_membership_n13_cdhit_0.90.tsv','\t')
ids={a for r in members if r['cluster_id'] in historical for a in r['members'].split(',')}
seqs={r['accession'].split('_domain')[0]:r['sequence'] for r in rows(ROOT/'reproducibility/data/refseq18_for_scoring.csv')}
actual={}
for a in ids:
    b=seqs[a]
    check(len(b)==len(q),'TG-1 comparison length '+a)
    diffs=','.join(f'{x}{i}{y}' for i,(x,y) in enumerate(zip(q,b),1) if x!=y) or '-'
    actual[a]=(sum(x==y for x,y in zip(q,b)),len(q),diffs)
archived=rows(t/'cluster4_identity.tsv','\t')
check(len(archived)==5 and actual=={r['accession']:(int(r['matches']),int(r['alignment_length']),r['substitutions']) for r in archived} and sorted(actual.values())==[(75,77,'V4I,A71V')]*2+[(77,77,'-')]*3,'TG-1: three exact matches, two 75/77 matches')
with tempfile.TemporaryDirectory() as d:
    out=Path(d)/'tg1.csv'
    subprocess.run([sys.executable,str(ROOT/'halc8_score_domains.py'),'--input',str(t/'query.csv'),'--output',str(out),'--input-type','domain','--ss-mode','none'],check=True,stdout=subprocess.DEVNULL)
    result=rows(out)
    check(out.read_bytes()==(t/'score.csv').read_bytes() and len(result)==1 and float(result[0]['total_score_9'])==6.0 and result[0]['verdict_v2']=='PROBABLE','TG-1 score: 6.0 PROBABLE')
print('Selected evidence checks: PASS (not full workflow reproduction)')
