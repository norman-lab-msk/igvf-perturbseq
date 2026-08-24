# RPE-1 Multiome Perturb-seq (CRISPRi) — analysis code

**Paper.** Metzner E*, Southard KM*, Norman TM (*equal contribution).
*Multiome Perturb-seq unlocks scalable discovery of integrated perturbation effects on the
transcriptome and epigenome.* Cell Systems **16**, 101161 (2025).
[10.1016/j.cels.2024.12.002](https://doi.org/10.1016/j.cels.2024.12.002)

## Screen

CRISPRi (dCas9-ZIM3) in hTERT RPE-1, gene expression and chromatin accessibility from the same
nuclei via 10x Multiome. 13 chromatin remodelers — SMARCE1, SMARCB1, ARID1A, SMARCA4, DPF2,
SMARCC1, SMARCC2, EP400, ACTL6A, DMAP1, SUZ12, EZH2, YY1 — plus 3 non-targeting controls, in one
channel (`Lane1_040`).

## The chain

    cellranger-arc output  +  guide assignment from the gRNA FASTQs
              |
              v
    matrices in which guide identity is a dimension
              |
              v
    element-level differential analyses

cellranger-arc has no CRISPR guide-capture modality, so the sgRNA library was sequenced separately
and called on its own. That is why the middle step exists: it is where the guide axis is
reconstructed and joined to the expression and accessibility data.

| Step | Directory | Produces |
|---|---|---|
| 0 | [`00_cellranger_arc/`](00_cellranger_arc/) | the joint matrix and the ATAC fragments — invocation and libraries file only, no lab code |
| 1 | [`01_guide_assignment/`](01_guide_assignment/) | per-(barcode, UMI, protospacer) table and the singlet calls |
| 2 | [`02_guide_assigned_matrices/`](02_guide_assigned_matrices/) | `rpe1_multiome_cell_by_guide.h5ad` — the guide axis — plus the singlet-filtered intermediates |
| 3 | [`03_differential_analyses/`](03_differential_analyses/) | the two element-level result tables, and the exact bytes submitted |

Step 1's code is **public in the paper repository**,
[`norman-lab-msk/multiomeperturbseq`](https://github.com/norman-lab-msk/multiomeperturbseq), and is
referenced rather than copied — see that step's README for the byte-level comparison against the
copy that was actually run. The normalisation and ATAC peak-calling that feed step 3 are in the same
repository, in `Fig1_vRevised.ipynb` and `Fig2_vRevised.ipynb`.

## "element" definition

The **perturbed promoter**, keyed on `intended_target_chr/_start/_end` — nearest EPD extended
promoter window to the guide. The two modalities differ only in readout: expression reports
`target_gene`, accessibility reports the peak's `chr/start/end`.

## Outputs

Run of record, 2026-08-21, in `results/`.

| File | Grain | Rows |
|---|---|---|
| `multiome_paper_guide_effect_matrix.csv` | element × gene | 14,868 |
| `multiome_paper_differential_peaks_by_guide.csv` | element × peak | 10,935 |

Columns — `effect_score`, `p_val`, `p_val_adj`, `guide_id`,
`intended_target_name`, `intended_target_chr/_start/_end` — plus `target_gene` (GEX) or
`chr/start/end` of the peak (ATAC).

`effect_score` is a different quantity in each: GEX is mean z-scored expression over the guide's
cells, in SD units; ATAC is log2 fold-change of mean accessibility against NTC, with a `1e-3`
pseudocount. See `03_differential_analyses/README.md` for why.

## Environment

`environment.yml` in this directory, on top of the repository baseline. **snapatac2 must be
2.6.0** — other versions do not reproduce the published Table S3. `perturbseq` is not packaged;
see [thomasmaxwellnorman/perturbseq_demo](https://github.com/thomasmaxwellnorman/perturbseq_demo).

## To run

Each step's README gives its own command. Paths default to the cluster locations used for the
published run; override them with the environment variables each step documents.
