# Step 3 — element-level differential analyses

Both modalities against the non-targeting controls, reported per **element** — the perturbed
promoter, the nearest EPD extended promoter window to the guide, keyed on
`intended_target_chr` / `_start` / `_end`.

## Contents

| File | Role |
|---|---|
| `multiome_data_organization.ipynb` | maps guides to promoter windows, then runs both differential analyses |
| `multiome_paper_igvf_guides.csv` | the 16-guide library table, with coordinates filled in by the notebook |
| `prepare_upload.py` | reformats the two result CSVs into the exact bytes submitted to IGVF |

Outputs land in `../results/`.

## The two analyses

**Expression.** `perturbseq.ks_de` on the normalised, regressed matrix: a Kolmogorov–Smirnov test of
each guide's cells against the NTC cells, Benjamini–Yekutieli control, effect score = the mean
z-scored expression of the measured gene in that guide's cells. Filtered to q < 0.1 →
14,868 (element, gene) calls.

**Accessibility.** Mann–Whitney U per guide against the NTC cells over the MACS3 peak matrix,
Benjamini–Hochberg control, effect score = log2 fold change of mean accessibility with a `1e-3`
pseudocount. Filtered to q < 0.1 → 10,935 (element, peak) calls.

Both are non-parametric tests of the same shape; the effect sizes differ because the matrices do.
Expression has already been normalised against the controls and regressed, so its values are
centred residuals on which a fold change is not defined; accessibility is sparse non-negative
paired-insertion counts, where a ratio of means is standard and a per-peak z-score would track
sparsity rather than the perturbation. **`effect_score` therefore means something different in each
output file** and the two are not comparable on a common scale.

## Reproducing the submitted files

```bash
jupyter nbconvert --execute multiome_data_organization.ipynb   # writes ../results/*.csv
python prepare_upload.py                                       # writes ../results/upload/*.tsv.gz
```

`prepare_upload.py` drops the pandas index column, rewrites coordinate columns as integers (asserted
lossless) and gzips with `mtime=0` so the md5 is reproducible. It changes no values. The run of
record: `…gex_results.tsv.gz` md5 `722a58d8a34a0cfa52d06cdf9a691be8`, `…atac_results.tsv.gz` md5
`0ba54097facd23e7361ad48caf86b881`.

## Environment

`../environment.yml`. **snapatac2 must be exactly 2.6.0** — other versions do not reproduce the
published ATAC results. `perturbseq` is not packaged; it is loaded from a local directory and is
pinned by a content digest rather than a release.
