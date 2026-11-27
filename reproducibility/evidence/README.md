# Supporting evidence for HalC8-Prioritizer v2.2

Release candidate — not yet approved for publication. This directory complements `reproducibility/verify.py`, which reproduces the three frozen scoring outputs. The scoring rules and profile HMM are not changed.

## Reproduction and scope

Run from the repository root:

```sh
python3 reproducibility/verify.py
python3 reproducibility/evidence/verify_evidence.py
```

The second command checks artifact hashes, the reconstructed RefSeq extraction (44 retained targets, 26 exact calibration matches, 18 remaining), three control scoring datasets and the corrected translated-CDS RNS21 scoring output. It also checks archived clustering counts and the locus inventory. These checks are not a new full HMMER, CD-HIT, MMseqs2, BLAST or Clinker run.

| Directory | Scope and limitations |
|---|---|
| calibration | Historical prior-label table and LOOCV/ablation scripts. Internal calibration, not external validation. Original scripts contain historical interpretations and approximate uncertainty calculations; consult the final manuscript for corrected interpretation and tie handling. |
| retrieval, extraction | Archived assembly inventory, CD-HIT and HMMER outputs; exact reconstructed domain extraction from the retained 44 source proteins. Full RefSeq snapshot and NR100 proteome are not bundled. Original HMM construction and extractor cannot be regenerated exactly from the supplied files. |
| controls | Three off-target proteome audits, 8 compatible hard negatives (plus ninth historical record), 2,062 genome-pilot records and HMMER checks of the four score-positive events. Original proteomes and pilot sampling not rerun. |
| comparators | Archived n39 mature-domain HMM, HMM-only ranking and precursor-model outputs, corrected translated-CDS RNS21 and Torres BLAST audit. Third-party models and Torres query FASTA are listed by hash in THIRD_PARTY_INPUTS.tsv, not distributed pending rights verification. |
| clustering | n18 and n13 archived memberships, n13 scripts and pairwise identity outputs; cluster ID mapping. CD-HIT/MMseqs/EMBOSS not rerun here. |
| sensitivity | Z3 and HG1 per-record evidence; Cluster 5 excluded from fully assessable shortlist without changing scorer. Incorrect historical AS7092 aggregate omitted; use locus-level inventory. |
| synteny | 54 historical records plus 2 post hoc records for the same C3 biological locus; 56 records (42 full, 14 truncated). Representative GenBank windows and archived Clinker links, not all historical genomic inputs. |
| figure_source | Corrected Figure 2 and 3 source TSVs/scripts; Figure 4 and S4/S5 sources are under synteny. Exact submitted figure correspondence and remaining Figure 1/S1-S3 provenance not verified. |

Historical files retain their original analytical names and relative paths. Some scripts require the original directory layout; use copies, not archived outputs, for reruns. Local user paths have been redacted where necessary. PROVENANCE.tsv identifies every selected original and distinguishes changed bytes from originals. Historical manifests may reference original hashes and paths. SHA256SUMS.txt describes this candidate, not the historical archive.

Historical NCBI Revision History screenshots underlying Table S10 were not recovered; this package does not independently verify their dated statuses. No experimental activity or unseen-population validation is claimed. TG-1 and unrelated exploratory projects are excluded pending confirmation of manuscript scope.

## Publication gate

Confirm exact submitted manuscript/supplement correspondence, remaining figure and genus-audit evidence, and both authors' permission before pushing. New materials are not automatically covered by the existing MIT or CC BY licence; verify attribution and redistribution rights, including the Chou-Fasman module and any third-party sources. If rights prevent distribution, document source, version, hash and retrieval or controlled access. Do not describe this archive as a complete end-to-end reproduction.
