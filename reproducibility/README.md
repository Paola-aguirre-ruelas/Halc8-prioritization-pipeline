# Core scoring reproducibility

This directory reproduces scoring on the fixed set of 18 calibration-set-novel RefSeq domains and 42 calibration-derived domains (no secondary structure and local Chou-Fasman configurations). It does not reconstruct the upstream proteome snapshot, the profile HMM, the historical full-mode calibration, or the other supplementary analyses. Score-positive labels are prioritization decisions, not experimentally measured activity.

From the repository root, install the requirements in the repository and run:

    python3 reproducibility/verify.py

The script verifies input and implementation hashes and compares all three output CSVs byte for byte against the archived outputs. It writes temporary outputs outside the repository. A difference in serialization or dependencies requires investigation; it does not automatically establish a scientific discrepancy.

The historical scorer SHA-256 was `0bd6dec9e75053f552c84585573be9b53671a00188bfcd7dadbfbd7f774f2e2b`. The public implementation has a different hash; its three output CSVs have been reproduced exactly against the archived tables using the same runtime and secondary-structure helper. Historical prior-label provenance: https://doi.org/10.5281/zenodo.21116990. The 42-sequence set is not independent validation.

Data provenance: the inputs and reference outputs are unchanged copies of the archived core scoring artifacts. See the data license and original biological database provenance. The original release scope does not explicitly cover the recovered Chou-Fasman helper; its redistribution license must be confirmed by the authors before public release.
