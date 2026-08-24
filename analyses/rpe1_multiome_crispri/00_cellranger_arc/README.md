# Step 0 — cellranger-arc

The joint GEX + ATAC count step. It is off-the-shelf 10x software, so there is no lab code here —
what is recorded is the exact invocation, taken from the run directory itself rather than from
anyone's notes.

## Invocation

`cellranger_arc_cmdline.txt`, verbatim from the run's `_cmdline`:

```
cellranger-arc count --id=Lane1_040 \
  --reference=/fscratch/eli/genomes/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
  --libraries=040_libraries.csv --localcores=40 --localmem=128
```

Versions, from the run's `_versions`: **cellranger-arc 2.0.2**, martian v4.0.5. The same version and
reference path are also printed in the header comments of `atac_fragments.tsv.gz`, so they can be
read back off the data without this directory.

## Libraries

`040_libraries.csv` is the file that invocation consumed: five FASTQ directories across two
modalities and five sequencing runs.

| Modality | Runs |
|---|---|
| Chromatin Accessibility | `DIANA_0696`, `DIANA_0699`, `DIANA_0700` |
| Gene Expression | `FAUCI_0130`, `FAUCI_0133` |

The T7 sgRNA library is deliberately absent: cellranger-arc has no CRISPR guide-capture modality, so
the guides were sequenced as their own library and called separately — that is step 1.

## Outputs used downstream

| Output | Used by | Submitted as |
|---|---|---|
| `filtered_feature_bc_matrix.h5` | steps 1, 2 | `cell by gene and peak matrix` (`IGVFFI6914SXML`) |
| `atac_fragments.tsv.gz` (+ `.tbi`) | step 2's peak calling | `fragments` |

Paths in both files are lilac paths, where the run was performed. The run directory has since been
copied to `/data1/collab005/202404_SIRLOIN_multiome_cellranger/Lane1_040/` on iris.
