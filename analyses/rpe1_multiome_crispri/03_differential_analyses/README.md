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
| `all_pairs_gex_results.py` | the expression analysis for every tested pair, with `log2_fold_change` and `significant` |

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

`all_pairs_gex_results.py` runs the expression test exactly as the notebook does and keeps all
13 × 9,992 = 129,896 (element, gene) pairs instead of the 14,868 at q < 0.1. Output:
`../results/rpe1_multiome_element_level_gex_results_all_pairs.tsv`, SHA-256
`f16634d83125a4a49769b32e2a8da278b1153b71c4360da4061085d287dbb70d`. It is the notebook's nine columns
in the submitted layout, then:

| Column | Meaning |
|---|---|
| `log2_fold_change` | `log2((CPM_guide + 1) / (CPM_NTC + 1))`, from raw counts summed over each guide's nuclei and over the 1,144 NTC nuclei in `040_singlets_gex_unprocessed.h5ad`. |
| `significant` | `TRUE` / `FALSE`, `p_val_adj < 0.1`. 14,868 rows are `TRUE`. |

Its `TRUE` rows, first nine columns, are the submitted expression table byte for byte.

`log2_fold_change` is descriptive. The p-values belong to the KS test behind `effect_score`, which
is taken after regression on mitochondrial fraction and cell-cycle scores; the fold change is not.
Among the `TRUE` rows the two agree in sign on 85%, and on 98% or more once `|effect_score| > 0.1`.
The fold change also shows what the test does not call: all 13 guides lower their own target
(log2 fold change −0.64 to −2.25), while only 6 of those 13 pairs reach q < 0.1.

26 rows — two gene symbols across the 13 elements — have an empty `target_gene` and no
`log2_fold_change`, because the notebook's symbol-to-Ensembl lookup does not resolve them. Five of
them are in the submitted table.

```bash
python all_pairs_gex_results.py --gex /writable/copy/of/gex_norm_regressed.hdf5
```

`perturbseq` opens its store read-write, so `--gex` must point at a copy you can write to; the file
is not modified (SHA-256 `2ba7f3ac585f8e2b04be423de8a8c2fd3a38a441c3f769a80061167f2865bd71` before and
after). Raw counts: SHA-256 `d13de8ff385f03c6cfe4188ca049ce2cb15075fabaa06b27f35a4e9d4d00db7a`.

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
