# HalC8-Prioritizer v2.2

First public release of HalC8-Prioritizer, implementing the frozen HalC8 v2.2 domain-scoring rules and including the corresponding custom profile HMM.

## Scope and interpretation

The public repository provides the frozen family-retrieval HMM and a standalone, rule-based domain scorer. The scorer accepts pre-extracted amino-acid domains; it does not itself reproduce the upstream proteome retrieval and HMM-based domain extraction workflow. Prioritization scores are not experimentally validated probabilities of antimicrobial activity. A score-positive verdict is a candidate-prioritization decision, not evidence of biological activity.

## Versioning

Version 2.2 is the first public release of HalC8-Prioritizer. The v2.2 identifier follows internal workflow-development numbering; earlier workflow iterations, including v1.x and v2.1, were not publicly released or separately registered. The bundled standalone scoring script retains its own component implementation identifier (`halc8_score_domains_v1.0`), which should not be interpreted as a separate public release of the complete workflow.

## Contents

- `halc8_score_domains.py` — standalone HalC8-family domain scoring script.
- `environment.yml` — conda environment specification.
- `requirements.txt` — Python package requirements.
- `LICENSE` — software license.
- `CITATION.cff` — citation metadata for this software repository.

## Installation

Using conda:

    conda env create -f environment.yml
    conda activate halc8-pipeline

Using pip:

    python3 -m pip install -r requirements.txt

## Basic usage

    python3 halc8_score_domains.py --input input_domains.csv --output scored_domains.csv --input-type domain --ss-mode none

Input CSV files require at least two columns: `accession` and `sequence`. With `--input-type domain`, the supplied sequence is scored directly. With `--ss-mode none`, secondary-structure scoring is disabled; results should not be interpreted as full-mode scores. The script writes a scoring CSV and an accompanying run manifest.

## License

See `LICENSE`.

## Data provenance

The 2021 prior-label table used for calibration-frame provenance is archived at Zenodo: https://doi.org/10.5281/zenodo.21116990. This table is used as a historical prior-label framework, not as an experimentally validated benchmark or independent external validation set.

## Frozen HMM profile

The frozen HalC8 v2.2 custom HMM profile used for family-level retrieval is included in `hmm/HalC8_custom_v2.2.hmm`; its SHA-256 checksum is provided in `hmm/HalC8_custom_v2.2.hmm.sha256`.

Verify the frozen HMM from the repository root:

    shasum -a 256 -c hmm/HalC8_custom_v2.2.hmm.sha256
