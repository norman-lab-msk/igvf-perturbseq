#!/usr/bin/env python3
"""Export the processed RPE-1 Multiome Perturb-seq expression matrix as AnnData.

The normalisation step (`Fig1_vRevised.ipynb` in `norman-lab-msk/multiomeperturbseq`) writes its
result through `perturbseq`'s `CellPopulation.to_hdf`, which is a pandas HDFStore with `matrix`,
`cell_list` and `gene_list` keys. That layout is readable only with the `perturbseq` library, which
is not packaged, so it is a poor thing to hand to anyone else.

This script rewrites exactly that object as a standard `.h5ad` — same values, same cells, same
genes, nothing recomputed:

    matrix     -> X        4,724 nuclei x 9,992 genes, normalised against non-targeting controls
                           within each GEM group and regressed on pct_counts_mito, S_score, G2M_score
    cell_list  -> obs      including `guide_identity`, the assigned perturbation
    gene_list  -> var      the gene filter's statistics

Float64 is narrowed to float32 on the way out. That is the one lossy step and it is deliberate: it
halves the file, and the values are z-scores whose seventh decimal place is noise. Pass
`--keep-float64` to disable it.
"""
import argparse
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).parent
SRC = Path(os.environ.get(
    "MULTIOME_GEX_NORMALIZED",
    "/data1/normantm/eli/T7/202404_SIRLOIN_multiome/share/gex_norm_regressed.hdf5"))


def export(src, keep_float64=False):
    with pd.HDFStore(src, "r") as store:
        matrix = store["matrix"]
        cells = store["cell_list"]
        genes = store["gene_list"]

    if not matrix.index.equals(cells.index):
        raise ValueError("matrix rows and cell_list are not aligned")
    if not matrix.columns.equals(genes.index):
        raise ValueError("matrix columns and gene_list are not aligned")

    X = matrix.to_numpy()
    if not keep_float64:
        X = X.astype(np.float32)

    adata = ad.AnnData(X=X, obs=cells.copy(), var=genes.copy())
    adata.uns["description"] = (
        "Normalised, regressed expression for the single-guide nuclei of RPE-1 Multiome Perturb-seq "
        "channel Lane1_040. Rewritten from the perturbseq CellPopulation HDFStore; values unchanged.")
    return adata


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--src", type=Path, default=SRC)
    ap.add_argument("--out", type=Path,
                    default=HERE / ".." / "results" / "rpe1_multiome_processed_gex.h5ad")
    ap.add_argument("--keep-float64", action="store_true")
    args = ap.parse_args()

    adata = export(args.src, args.keep_float64)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out, compression="gzip")

    print(f"cells {adata.n_obs}  genes {adata.n_vars}  dtype {adata.X.dtype}")
    print(f"  obs columns: {', '.join(adata.obs.columns)}")
    print(f"  perturbations: {sorted(adata.obs['guide_identity'].unique())}")
    print(f"  wrote {args.out} ({args.out.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
