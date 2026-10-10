# Hs27 fibroblast CRISPRa Perturb-seq — element-level results

Converts the screen's mean-population matrix into element-level table. Results of record:
2026-03-03.

**Paper.** Southard KM*, Ardy RC*, Tang A, O'Sullivan DD, Metzner E, Guruvayurappan K, Norman TM. (*equal contribution)
*Comprehensive transcription factor perturbations recapitulate fibroblast transcriptional states.*
Nature Genetics **57**, 2323–2334 (2025).
[10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1)

## "element" definition

A **promoter window**, keyed on
`(intended_target_name, intended_target_chr, intended_target_start, intended_target_end)`. A gene
with two annotated windows stays two elements — 307 elements map onto 291 genes.

Each element is represented by a **single guide**, the one with the highest `de_genes_fibro` among
its hit guides.

## Upstream

The input matrix comes from
**[norman-lab-msk/TFs_CRISPRa](https://github.com/norman-lab-msk/TFs_CRISPRa)** — Cell Ranger 7.1.0
on GRCh38-2020-A, guides at a 5-UMI threshold, per-gem-group normalization against control cells,
occupancy corrected regression on the guide design matrix.

The result is `fibroblast_CRISPRa_mean_pop.h5ad` (10,916 guides × 4,914 genes): `X` holds effect
scores; layers `p`, `adj_p`, `masked` hold p-values, adjusted p-values and target-masked
coefficients.

## Processing

`generate_element_level_results.py`:

1. Joins Table S4 (`Guide Activity`) onto `adata.obs` by `guide_identity`.
2. Selects hit guides — `seed_driven_fibro == False AND bad_seed == False AND (active_fibro OR
   expanded_active_fibro)`.
3. Picks one representative guide per element — highest `de_genes_fibro`.
4. Melts effect scores and both p-value layers to one row per element × gene, maps symbols to
   unversioned Ensembl IDs, filters to `p_val_adj < 0.05`.

`add_log2_fold_change.py` then extends the same table to **every tested pair** — step 4 without the
filter — and adds two columns:

- `log2_fold_change`, from the count-scale regression of the same screen (TFs_CRISPRa, *Step 4 -
  Simplified regressions*): raw UMI counts regressed on the guide design matrix, so the `control`
  coefficient is baseline UMIs per cell and each guide coefficient is the change from it. Modelled
  expression with the guide is `control + coefficient`, clipped at zero, and the column is
  `log2((CPM + 1) / (CPM_control + 1))`.
- `significant`, `TRUE` where `p_val_adj < 0.05` — the rows of `element_level_results.tsv`.

## Inputs

Pinned by SHA-256, paths relative to this directory.

| File | Shipped | SHA-256 |
|---|---|---|
| `merged_guides_promoters_UPSTREAM_PRIORITIZED.csv` | yes | `1ac915099d18032148841b862331ce7bb33d24bf811f9bee01424c4b534c6f0f` |
| `Supplementary_Table1.xlsx` — 17 MB, distributed with the paper. Only sheet `Table S4 Guide Activity` is read. | no | `dc125f4727d650163fbfb08175d2a48ea77406787e77d3fa21a134639266a6c0` |
| `fibroblast_CRISPRa_mean_pop.h5ad` — 1.7 GB, from Zenodo below. | no | `6e10a2605a3dbc448ac4756399adc947a493a75f8961d3793af8f851aa0411cf` |
| `20240310_fibroblast_final_low_UMI_count_regressions_coef_mean_pop.h5ad` — 12 GB, 10,917 × 30,395. `X` is the coefficient matrix TFs_CRISPRa Step 4 writes as `…coef_mean_pop.hdf`, converted to AnnData. Read only by `add_log2_fold_change.py`. | no | `000e4a776d016b66e8da5ce5d0aebf137406e016fc2d83683e518455b197e2a4` |

**Data availability.** Mean-population matrix:
[10.5281/zenodo.15200179](https://doi.org/10.5281/zenodo.15200179) (byte-identical to the run of
record, md5 `ba44c7813903bb5df900348d6b0d589a`). Other processed data:
<https://zenodo.org/communities/normanlabmsk/>.

The count-scale regression matrix is not deposited. It is the output of *Step 4 - Simplified
regressions for newly expressed genes and on-target activation in Hs27 experiment* in
[norman-lab-msk/TFs_CRISPRa](https://github.com/norman-lab-msk/TFs_CRISPRa), run on two files from
the same Zenodo record: `fibroblast_CRISPRa_final_pop.h5ad` (md5 `7abb92c81896ec3334d144a99dd544d1`)
and `fibroblast_CRISPRa_aggr_total_guide_umis.h5` (md5 `3d6f1c7cf3ac5a3c034d67ca6525e587`). The
notebook writes `…coef_mean_pop.hdf`; its `/matrix` table is `X` in the AnnData read here. Raw FASTQs: SRA
[PRJNA1108254](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1108254).

## Outputs

`results/element_level_results.tsv` — one row per element × gene at `p_val_adj < 0.05`.

| Column | Meaning |
|---|---|
| `effect_score` | Deviation from control in SD units. Negative is down. |
| `p_val` / `p_val_adj` | Regression p-value, and Benjamini–Hochberg adjusted. |
| `guide_id` | Representative guide for the element, `SYMBOL_PROTOSPACER`. |
| `target_gene` | Ensembl ID of the **measured** gene whose expression changed. |
| `intended_target_name` | Ensembl ID of the gene whose promoter was **perturbed**. |
| `intended_target_chr/_start/_end` | Promoter window, always 1,500 bp wide. UCSC-style chromosome (`chr4`). |

`results/element_level_results_guide_summary.tsv` — one row per element, with `n_de_genes`,
`n_upregulated`, `n_downregulated` and the element coordinates.

`results/element_level_results_all_pairs.tsv` — one row per tested element × gene: 1,538,082 rows,
313 elements × 4,914 genes. 257 MB, so it is not committed; SHA-256
`d27c1b19ee34c82e72369fef5512967ecc384adab3a4b3e619dc8de700238e8c`. The columns above, then:

| Column | Meaning |
|---|---|
| `log2_fold_change` | Model-based log2 fold change against control, 1 CPM pseudocount. Descriptive: the p-values belong to `effect_score`. |
| `significant` | `TRUE` / `FALSE`, `p_val_adj < 0.05`. 31,047 rows are `TRUE`. |

Its `TRUE` rows, first nine columns, are `element_level_results.tsv` byte for byte. It has six more
elements than that file: active elements with no gene at `p_val_adj < 0.05`.

`log2_fold_change` is a ratio of regression estimates, not of two observed groups of cells — the
droplets are overloaded, so no cell carries one guide alone. Where the additive model predicts
expression at or below zero (14,768 rows, 1.0%; 317 of them `TRUE`) the value is
`log2(1 / (CPM_control + 1))`: it says the gene is modelled as off and scales with the gene's
baseline, not with the strength of repression. These sit on low-baseline genes.

## To run

Fetch the matrix from Zenodo, then run **from this directory**

```bash
python generate_element_level_results.py --out results/element_level_results.tsv

# every tested pair, with log2_fold_change and significant
python generate_element_level_results.py --padj_thresh 2 --out element_level_results_unfiltered.tsv
python add_log2_fold_change.py
```
Tested with Python 3.13.7, scanpy 1.11.4, anndata 0.12.2, pandas 2.3.2, numpy 2.3.3, scipy 1.16.2,
h5py 3.14.0, openpyxl 3.1.5.
