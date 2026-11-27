#!/usr/bin/env python3
"""Render Figure 4 from archived Major7 ZIP. Coordinates/strands from GenBank;
colour assignments from direct archived Clinker links. Table S7 contextual
interpretation. No synthetic CDS. Display orientation normalized per candidate.
Usage: python generate_figure4_synteny.py HalC8_Major7_evidence.zip [output_dir]
Dependencies: biopython, matplotlib.
"""
import sys, zipfile, json, csv, io, hashlib
from pathlib import Path
from itertools import combinations
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch

archive=zipfile.ZipFile(sys.argv[1]); out=Path(sys.argv[2]) if len(sys.argv)>2 else Path('.')
out.mkdir(parents=True,exist_ok=True)
root='synteny_figures_final/'
refs=['halU','halR','proC8','halT1','halT2','halT3']
panels={}
for panel in ['S4_AS7092_supported_C1_C3','S5_non_AS7092_C2_C4_C5']:
    path=root+panel+'/'
    rows=list(csv.DictReader(io.StringIO(archive.read(path+'manifest.tsv').decode()),delimiter='\t'))
    data=json.loads(archive.read(path+'clinker_outputs/clinker.json'))['data']
    tracks={}
    for row in rows:
        b=archive.read(path+'inputs/'+row['gbk_filename']); assert hashlib.sha256(b).hexdigest()==row['sha256']
        record=SeqIO.read(io.StringIO(b.decode()),'genbank')
        genes=[f for f in record.features if f.type=='CDS']
        acc=row['candidate_accession']; c=next(c for c in data['clusters'] if c['name']==row['display_name'])
        cg=c['loci'][0]['genes']
        target=next(f for f in genes if (f.qualifiers.get('gene')==['proC8'] if acc=='-' else acc in f.qualifiers.get('protein_id',[])))
        tracks[acc]=(row,record,genes,cg,target)
    links={frozenset((x['query']['uid'],x['target']['uid'])):x['identity'] for x in data['links']}
    ref={g['uid']:g['names']['gene'] for g in tracks['-'][3]}
    panels[panel[:2]]=tracks,links,ref

# Derive C4 core from adjacent CDS in Clinker orientation, then verify all six
# pairwise links for each of the five neighbors and the candidate itself.
tracks,links,_=panels['S5']; core={}
for acc in ['WP_006600928.1','WP_004043226.1','WP_115805229.1','WP_004975848.1']:
    gs=sorted(tracks[acc][3],key=lambda g:g['start']); i=next(i for i,g in enumerate(gs) if g['label']==acc)
    core[acc]=gs[i-1:i+5]
for i in range(6):
    assert all(links.get(frozenset((a[i]['uid'],b[i]['uid'])),0)>=.3 for a,b in combinations(core.values(),2))
colors={'candidate':('#EF9F27','#854F0B'),'linked':('#5DCAA5','#0F6E56'),'neighbor':('#F0997B','#993C1D'),'other':('#B4B2A9','#5F5E5A')}
fig,axes=plt.subplots(2,1,figsize=(9.5,7.1),gridspec_kw={'height_ratios':[3.4,3.5]})
fig.subplots_adjust(left=.07,right=.96,top=.95,bottom=.17,hspace=.35)
selections=[('S4',['-','WP_238705315.1','WP_015299281.1'],'A   AS7092-supported cassette context'),('S5',['WP_006600928.1','WP_004043226.1','WP_136361253.1'],'B   Score-positive, non-AS7092 contexts')]
for ax,(p,accs,title) in zip(axes,selections):
    tracks,links,ref=panels[p]; ax.set_title(title,loc='left',fontsize=13,fontweight='bold',pad=10)
    ax.set_xlim(-5700,5700); ax.set_ylim(-.05,3.05); ax.set_yticks([]); ax.set_xticks([-5000,0,5000],['−5 kb','0','+5 kb']); ax.tick_params(axis='x',labelsize=8,length=3)
    for side in ['top','right','left']: ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_color('#AAAAAA'); ax.axvline(0,color='#BBBBBB',linestyle=':',linewidth=.8,zorder=0)
    for idx,acc in enumerate(accs):
        row,rec,genes,cg,target=tracks[acc]; y=2.4-idx*.98
        anchor=(int(target.location.start)+int(target.location.end))/2; orient=target.location.strand
        label='AS7092 reference · DQ137141.1' if acc=='-' else acc+'   '+row['display_name'].split(' . ',2)[2]+'   · '+row['cluster']
        ax.text(-5600,y+.3,label,fontsize=9.2,va='bottom')
        linked={g['label']:next((name for uid,name in ref.items() if links.get(frozenset((uid,g['uid'])),0)>=.3),None) for g in cg}
        recurrent={g['label'] for i,g in enumerate(core.get(acc,[])) if i!=1}
        for f in genes:
            a,b=sorted([(int(f.location.start)-anchor)*orient,(int(f.location.end)-anchor)*orient]); strand=f.location.strand*orient
            if a < -5700 or b > 5700: continue # only full depicted CDS within display frame
            pid=f.qualifiers.get('protein_id',f.qualifiers.get('locus_tag',['']))[0]
            gene=f.qualifiers.get('gene',[''])[0]; lname=linked.get(pid)
            cat='candidate' if f is target else 'linked' if ((acc=='-' and gene in refs) or (p=='S4' and lname and lname!='proC8')) else 'neighbor' if pid in recurrent else 'other'
            fc,ec=colors[cat]; tip=min(120,(b-a)*.4); h=.14
            if strand==1: pts=[(a,y-h),(b-tip,y-h),(b,y),(b-tip,y+h),(a,y+h)]
            else: pts=[(b,y-h),(a+tip,y-h),(a,y),(a+tip,y+h),(b,y+h)]
            ax.add_patch(Polygon(pts,closed=True,facecolor=fc,edgecolor=ec,linewidth=.7))
            if acc=='-' and gene in refs: ax.text((a+b)/2,y-.19,gene+('*' if gene=='halT3' else ''),fontsize=8,ha='center',va='top')
        # Labels refer to the actual candidate-centered coordinate frame.
    if p=='S5': ax.text(-5600,-.02,'',fontsize=8)
fig.legend(handles=[Patch(facecolor=colors[k][0],edgecolor=colors[k][1],label=l) for k,l in [('candidate','proC8 / candidate proC8-like CDS'),('linked','AS7092-linked cassette CDS'),('neighbor','Conserved Cluster 4 neighbouring CDS'),('other','Other annotated CDS')]],loc='lower center',bbox_to_anchor=(.5,.068),ncol=2,frameon=False,fontsize=9)
fig.text(.07,.045,'Candidate-centred coordinates; tracks oriented with candidate CDS forward.\nOnly CDS within the displayed interval are drawn. * Reference halT3 is partial at the record boundary.',fontsize=8,color='#444444')
for ext in ['pdf','svg','png']: fig.savefig(out/('Figure_4_verified.'+ext),dpi=300,facecolor='white')
