# Step 2 — guide-assigned matrices

Turns the cellranger-arc output plus the step 1 guide calls into the objects the differential
analyses run on. This is the step where **guide identity stops being a label and becomes a
dimension**.

## Contents

| File | Produces | Fate |
|---|---|---|
| `build_guide_matrix.py` | `rpe1_multiome_cell_by_guide.h5ad` — 9,318 nuclei × 16 guides, sgRNA UMI counts | submitted to IGVF as `cell by guide matrix` |
| `singlet_packaging.ipynb` | `040_singlets_gex_unprocessed.h5ad` and `singlets_assigned.csv` | intermediates — not submitted, but `singlets_assigned.csv` is where `build_guide_matrix.py` gets `guide_target`, `called` and `n_guides` |

`040_singlets_gex_unprocessed.h5ad` is deliberately not submitted: its matrix is a strict subset of
the cellranger-arc matrix already on the portal (the called barcodes × the gene features), so it
would be a second copy of data that is already there.

`singlet_packaging.ipynb` is Eli's `data_request.ipynb` with outputs stripped. It filters to the
4,724 nuclei carrying exactly one called guide and writes the raw counts for them, with the call and
the joint GEX/ATAC quality metrics in `obs`. **It has not been trimmed to just the steps that
produce those two files** — do that with Eli before relying on it as a one-command reproduction.

## Running

```bash
python build_guide_matrix.py                 # writes ../results/rpe1_multiome_cell_by_guide.h5ad
python build_guide_matrix.py --out /path/to/file.h5ad
```

Inputs default to the cluster locations of the published run and can be overridden with
`MULTIOME_GUIDE_UMIS`, `MULTIOME_CELLRANGER_H5`, `MULTIOME_SINGLET_CALLS` and
`MULTIOME_GUIDE_LIBRARY`. The UMI table lives on lilac; everything else is on iris.

The build is deterministic: the run of record produces md5 `8d97cbd46e76ce61df0b8f507d96ac1c`,
613,937 bytes, 67,334 non-zero entries.

## Two things worth knowing

**The three non-targeting guides are separate columns here.** Everything downstream collapses them
to a single `NTC` label, so this matrix is the only place in the submission where NTC identity is
resolved to the individual guide.

**`n_guides_detected` is not `n_guides`.** The first counts guides with at least one UMI in the
nucleus — raw evidence. The second is the caller's decision after thresholding. They disagree often,
which is the point of keeping the counts rather than only the call.
