# Retrospective TG-1 comparison

Source: Xing et al., Protein Expression and Purification 242 (2026), 106954, https://doi.org/10.1016/j.pep.2026.106954 (Figure 2; sequence transcribed in the archived 2026-11-26 analysis). The query is identical to three already distributed RefSeq mature domains. The publication itself is not redistributed here. The transcription has not been independently rechecked against the original figure in this verification.

`query.csv` and `score.csv` are unchanged archived inputs and outputs. `cluster4_identity.tsv` is the five-record subset of the archived identity analysis, expressed as ungapped residue comparisons with substitutions relative to TG-1. The verifier checks all five against the existing public RefSeq18 inputs and cluster mapping, and reruns the frozen scorer with `--input-type domain --ss-mode none`. The result is 6.0 / 9.0, PROBABLE. Three domains match 77/77; two match 75/77 (V4I and A71V).

The archived HMMER search predates the cited publication; the final scoring manifest postdates it. This is retrospective sequence concordance, not prospective or independent predictive validation, nor evidence that every matching candidate has experimentally verified activity. The broader historical HMM30 identity analysis is outside this verification.

## Source integrity

Original private analysis artifacts (SHA-256):
- `query.csv`: `84adde56581069fcaaed6b55abb01bb12f98394116a536a45a32bc4dccf8a4bf`
- `score.csv`: `b075f9d9e83724a39a0e445774fcbed99604caed204a3a57df3e86fab80724a2`
- `identity.tsv`: `6ccc6b29d832288ccc198ba6034927b00e37c61d4cbc414b55db95bd09d3363d`
- `score.manifest.json`: `d884f9756e0608ba24c9a931e772ada8e56ca89b93889eefb255bad9c2488ec5`
- `tier1_archived_rows.txt`: `6f40d71fc22a0dc8506e92e1b866cd543cafa28dd9ab8c0a47634aee34de107f`

Original scoring run recorded 2026-11-26T19:31:07Z; archived scorer SHA-256 `0bd6dec9e75053f552c84585573be9b53671a00188bfcd7dadbfbd7f774f2e2b`. The equivalent public scorer is used for the reproduction.
