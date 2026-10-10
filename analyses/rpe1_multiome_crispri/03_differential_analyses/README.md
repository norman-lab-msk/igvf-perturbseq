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
| `all_pairs_gex_results.py` | the expression analysis for every tested pair, with `log2_fold_change` |

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

## Every tested pair

`all_pairs_gex_results.py` runs the notebook's expression test without the q < 0.1 filter:
13 × 9,992 = 129,896 pairs, written to
`../results/rpe1_multiome_element_level_gex_results_all_pairs.tsv`. Its rows at `p_val_adj < 0.1`
are the submitted table. Same columns, plus:

| Column | Meaning |
|---|---|
| `log2_fold_change` | `log2((CPM_guide + 1) / (CPM_NTC + 1))` from raw counts summed over each guide's nuclei and the NTC nuclei. |

```bash
python all_pairs_gex_results.py --gex /writable/copy/of/gex_norm_regressed.hdf5
```

`perturbseq` opens its store read-write, so `--gex` needs a writable copy.

## Reproducing the submitted files

```bash
jupyter nbconvert --execute multiome_data_organization.ipynb   # writes ../results/*.csv
python prepare_upload.py                                       # writes ../results/upload/*.tsv
```

`prepare_upload.py` drops the pandas index column and rewrites coordinate columns as integers
(asserted lossless). It changes no values. The files were handed over uncompressed: `…gex_results.tsv`
md5 `219be77954867da29c48deb8178278d1`, `…atac_results.tsv` md5 `c912e674dad64fe9de18c8b8bbf132ce`.
The portal stores them gzipped and records these as `content_md5sum`; the `md5sum` it reports
(`de0e31a4…`, `60df2f12…`) is of its own compressed copy. `--gzip` writes a compressed form with
`mtime=0`, so that md5 is reproducible too.

## Environment

`../environment.yml`. **snapatac2 must be exactly 2.6.0** — other versions do not reproduce the
published ATAC results. `perturbseq` is not packaged; it is loaded from a local directory and is
pinned by a content digest rather than a release.
