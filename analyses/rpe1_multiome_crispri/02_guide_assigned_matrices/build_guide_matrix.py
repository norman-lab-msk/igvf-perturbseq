#!/usr/bin/env python3
"""Build the cell-by-guide matrix for the RPE-1 Multiome Perturb-seq channel Lane1_040.

Step 2 of the chain: cellranger-arc output + guide calls -> a matrix in which guide identity is a
dimension. cellranger-arc has no CRISPR guide-capture modality, so the T7 sgRNA library was
sequenced and called on its own (step 1) and guide identity survives only as a per-cell label. This
script recovers it as an axis, from the guide caller's own per-UMI output.

Inputs
------
UMI table   `Lane1_040_filtered_barcode_umis.csv.gz` from step 1 — one row per (cell barcode, UMI,
            protospacer). `identity` is at full spacer resolution, so the three non-targeting guides
            are distinct here even though the call collapses them to `NTC`.
cellranger  `filtered_feature_bc_matrix.h5`, for the authoritative cell-barcode set and order.
calls       `singlets_assigned.csv` from step 2's packaging notebook, for the per-nucleus call.
library     `multiome_paper_igvf_guides.csv` (step 3), for each guide's intended target.

Output
------
`rpe1_multiome_cell_by_guide.h5ad`
  X    cells x guides, uint32, distinct sgRNA UMIs supporting each guide in each cell
  obs  every cellranger-arc filtered barcode, with `guide_target`, `called`, `n_guides`,
       `total_guide_umis`, `n_guides_detected`
  var  one row per guide: `guide_id`, `spacer`, `targeting`, `intended_target_name`

Paths default to the cluster locations used for the published run; override with the environment
variables named below, or with --out for the output.

Barcode note: the UMI table has the `-1` suffix stripped; it is restored here so barcodes match the
expression and accessibility matrices exactly.
"""
import argparse
import os
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

HERE = Path(__file__).parent

UMIS = Path(os.environ.get(
    "MULTIOME_GUIDE_UMIS",
    "/data/norman/eli/T7/202404_SIRLOIN_multiome/guide_calling/Lane1_040_filtered_barcode_umis.csv.gz"))
CELLRANGER = Path(os.environ.get(
    "MULTIOME_CELLRANGER_H5",
    "/data1/collab005/202404_SIRLOIN_multiome_cellranger/Lane1_040/outs/filtered_feature_bc_matrix.h5"))
CALLS = Path(os.environ.get(
    "MULTIOME_SINGLET_CALLS",
    "/data1/normantm/eli/T7/202404_SIRLOIN_multiome/share/singlets_assigned.csv"))
GUIDES = Path(os.environ.get(
    "MULTIOME_GUIDE_LIBRARY",
    HERE / ".." / "03_differential_analyses" / "multiome_paper_igvf_guides.csv"))


def build(umis_path, cellranger_path, calls_path, guides_path):
    umis = pd.read_csv(umis_path)
    umis["CB"] = umis["CB"].astype(str) + "-1"

    with h5py.File(cellranger_path, "r") as h:
        barcodes = [b.decode() for b in h["matrix"]["barcodes"][:]]

    counts = (umis.groupby(["CB", "identity"], observed=True)["UB"]
                  .nunique()                       # UMIs, not reads
                  .rename("umis").reset_index())

    guide_ids = sorted(counts["identity"].unique())
    row = {b: i for i, b in enumerate(barcodes)}
    col = {g: i for i, g in enumerate(guide_ids)}

    keep = counts["CB"].isin(row)
    dropped = int((~keep).sum())
    counts = counts[keep]

    X = csr_matrix(
        (counts["umis"].to_numpy(dtype=np.uint32),
         (counts["CB"].map(row).to_numpy(), counts["identity"].map(col).to_numpy())),
        shape=(len(barcodes), len(guide_ids)), dtype=np.uint32)

    obs = pd.DataFrame(index=pd.Index(barcodes, name=None))
    calls = pd.read_csv(calls_path, index_col=0)
    obs["guide_target"] = calls["guide_target"].reindex(obs.index).fillna("").astype(str)
    obs["called"] = calls["called"].reindex(obs.index).fillna(False).astype(bool)
    obs["n_guides"] = calls["n_guides"].reindex(obs.index).fillna(0).astype(int)
    obs["total_guide_umis"] = np.asarray(X.sum(axis=1)).ravel().astype(int)
    obs["n_guides_detected"] = np.diff(X.indptr).astype(int)

    by_id = pd.read_csv(guides_path, index_col=0).set_index("guide_id")
    var = pd.DataFrame(index=pd.Index(guide_ids, name=None))
    var["guide_id"] = guide_ids
    var["spacer"] = [g.split("_")[-1] for g in guide_ids]
    var["targeting"] = [not g.startswith("NTC") for g in guide_ids]
    var["intended_target_name"] = [
        by_id["intended_target_name"].get(g, "") if g in by_id.index else "" for g in guide_ids]

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.uns["description"] = (
        "Per-cell guide UMI counts for RPE-1 Multiome Perturb-seq channel Lane1_040, from the T7 "
        "sgRNA library. Guide identity is the var axis; the downstream singlet call is obs.guide_target.")
    return adata, dropped


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=HERE / ".." / "results" / "rpe1_multiome_cell_by_guide.h5ad")
    ap.add_argument("--umis", type=Path, default=UMIS)
    ap.add_argument("--cellranger", type=Path, default=CELLRANGER)
    ap.add_argument("--calls", type=Path, default=CALLS)
    ap.add_argument("--guides", type=Path, default=GUIDES)
    args = ap.parse_args()

    adata, dropped = build(args.umis, args.cellranger, args.calls, args.guides)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out, compression="gzip")

    obs = adata.obs
    print(f"cells {adata.n_obs}  guides {adata.n_vars}  nonzero {adata.X.nnz}")
    print(f"  cells with >=1 guide UMI : {int((obs['total_guide_umis'] > 0).sum())}")
    print(f"  cells with a singlet call: {int(obs['called'].sum())}")
    print(f"  UMI rows on barcodes outside the filtered set, dropped: {dropped}")
    print(f"  wrote {args.out} ({args.out.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
