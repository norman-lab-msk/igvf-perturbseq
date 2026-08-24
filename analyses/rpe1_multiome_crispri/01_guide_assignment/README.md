# Guide assignment

The step that turns the T7 sgRNA library into a per-nucleus guide call. It is upstream of
`../notebooks/multiome_data_organization.ipynb`, which starts from objects that already carry
`guide_identity`.

**The code is not copied here — it is already public**, in the paper repository
[`norman-lab-msk/multiomeperturbseq`](https://github.com/norman-lab-msk/multiomeperturbseq):
`cropseq_call_guides_multiome.py` and its dependency `isotonic_barcode_calling.py`. Pin to a commit;
that repository has no tags.

Two notes for anyone reproducing the published run:

* `isotonic_barcode_calling.py` in that repository is byte-identical (md5 `5f5d81db…`) to the copy
  that was run from `/data1/normantm/eli/software`.
* `cropseq_call_guides_multiome.py` differs from the copy that was run by 15 lines: a `chdir`
  save/restore and hexbin plot styling. No calling logic differs.

## What it does

1. Takes the GEX-filtered cell barcodes from the cellranger-arc output.
2. Runs cellranger on the T7 guide FASTQs for a cell-barcode-labelled BAM.
3. Extracts the protospacer per read and collapses reads to (cell barcode, UMI, protospacer),
   writing `{sample}_raw_barcode_umis.csv.gz` and `{sample}_filtered_barcode_umis.csv.gz`.
4. Calls one guide per nucleus from those UMI counts by isotonic regression, keeping nuclei with
   exactly one confident call.

For `Lane1_040` the filtered UMI table is 314,922 (barcode, UMI, protospacer) rows over 9,251 of the
9,318 filtered nuclei, and the call retains **4,724** single-guide nuclei. `identity` is at full
spacer resolution there, so the three non-targeting guides stay distinct; only the call collapses
them to `NTC`.

## Output consumed downstream

`Lane1_040_filtered_barcode_umis.csv.gz` — the input to
[`../02_guide_assigned_matrices/build_guide_matrix.py`](../02_guide_assigned_matrices/), which turns
it into the matrix where guide identity is a dimension.
