#!/usr/bin/env python3
"""Expression results for every tested (element, gene) pair of the RPE-1 Multiome Perturb-seq screen.

`multiome_data_organization.ipynb` runs the expression test and keeps only the pairs at q < 0.1.
This script runs the same test on the same matrix and keeps all of them, so a pair that was tested
and not called is visible too. It adds one column to the notebook's nine:

  log2_fold_change   log2((CPM_guide + 1) / (CPM_NTC + 1)) from the raw counts of the same nuclei,
                     summed over each guide's nuclei and over the NTC nuclei. Descriptive only:
                     `p_val` and `p_val_adj` belong to the test behind `effect_score`, and the two
                     can disagree in sign where the effect is small.

The test itself is unchanged: `perturbseq.ks_de` of each guide's cells against the NTC cells on the
normalised, regressed matrix, Benjamini-Yekutieli within each guide, `effect_score` = mean z-scored
expression. Output is in the submitted layout (tab-separated, no index column, integer
coordinates); its rows at p_val_adj < 0.1 are the notebook's rows.

`target_gene` differs from the notebook's for seven genes. The notebook looks the Ensembl ID up by
gene name, which leaves two de-duplicated symbols blank and gives three genes that share a name with
another the other gene's ID; here the ID comes from Cell Ranger's feature table. Two more carry an
ID that is not in GENCODE v43 and are given their v43 ID (`../../gencode_v43_gene_ids.tsv`).
"""
import argparse
import os
import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).parent
DATA_ROOT = Path(os.environ.get("MULTIOME_DATA_ROOT", "/data1/normantm/eli"))
PERTURBSEQ = Path(os.environ.get("PERTURBSEQ_PATH", DATA_ROOT / "software"))
RESULTS = Path(os.environ.get("MULTIOME_RESULTS", HERE / ".." / "results")).resolve()

SHARE = DATA_ROOT / "T7" / "202404_SIRLOIN_multiome" / "share"
GEX_HDF5 = SHARE / "gex_norm_regressed.hdf5"
RAW_H5AD = SHARE / "040_singlets_gex_unprocessed.h5ad"
ID_UPDATES = HERE / ".." / ".." / "gencode_v43_gene_ids.tsv"
GUIDES_CSV = HERE / "multiome_paper_igvf_guides.csv"
OUT = RESULTS / "rpe1_multiome_element_level_gex_results_all_pairs.tsv"

Q_THRESH = 0.1
PSEUDOCOUNT = 1.0   # CPM
COLUMNS = ["effect_score", "p_val", "p_val_adj", "guide_id", "target_gene", "intended_target_name",
           "intended_target_chr", "intended_target_start", "intended_target_end"]


def melt(frame, value_name):
    long = (frame.reset_index().rename({"index": "gene"}, axis=1)
            .melt(id_vars="gene", var_name="guide_target", value_name=value_name))
    long.index = long.gene + "_" + long.guide_target
    return long


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    # perturbseq opens the store read-write, so point this at a writable copy if the original is not.
    p.add_argument("--gex", type=Path, default=GEX_HDF5)
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--n_jobs", type=int, default=12)
    args = p.parse_args()

    sys.path.append(str(PERTURBSEQ))
    from perturbseq import CellPopulation, ks_de

    guides = pd.read_csv(GUIDES_CSV, index_col=0)
    elements = guides[["guide_id", "gene_symbol", "intended_target_name", "intended_target_chr",
                       "intended_target_start", "intended_target_end"]]

    pop = CellPopulation.from_hdf(str(args.gex))
    pop.genes["gene_name"] = pop.genes.index
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        meanpop = pop.average("guide_identity")
        ks, ps, qs = ks_de(pop, key="guide_identity", control_cells="guide_identity == 'NTC'",
                           n_jobs=args.n_jobs)

    tested = melt(meanpop.matrix.T, "zscore")
    tested = tested.join(melt(ps.drop("NTC", axis=1), "p").p).join(melt(qs.drop("NTC", axis=1), "q").q)
    tested = tested.dropna(subset=["p"])          # the NTC column, which is not tested against itself

    res = (tested.reset_index(drop=True)
           .rename({"gene": "target_gene", "zscore": "effect_score", "p": "p_val", "q": "p_val_adj",
                    "guide_target": "gene_symbol"}, axis=1)
           .merge(elements, on="gene_symbol"))
    # Readout genes are named by Cell Ranger's symbols, which are de-duplicated (`MATR3-1`), so the
    # Ensembl ID is read from Cell Ranger's own feature table rather than looked up by name.
    raw = ad.read_h5ad(RAW_H5AD)
    cellranger_id = raw.var["gene_ids"]
    res["target_gene"] = res["target_gene"].map(cellranger_id)
    assert res["target_gene"].notna().all()
    res = res[COLUMNS].copy()
    for column in ("intended_target_start", "intended_target_end"):
        assert (res[column] % 1 == 0).all(), f"{column} is not integral"
        res[column] = res[column].astype(int)

    group = raw.obs["guide_target"].astype(str).values
    cpm = {}
    for name in np.unique(group):
        counts = np.asarray(raw.X[group == name].sum(axis=0)).ravel()
        cpm[name] = counts / counts.sum() * 1e6
    gene = pd.Index(raw.var["gene_ids"]).get_indexer(res["target_gene"])
    target = res["guide_id"].str.split("_").str[0].values
    log2fc = np.full(len(res), np.nan)
    for name in np.unique(target):
        rows = (target == name) & (gene >= 0)
        log2fc[rows] = np.log2((cpm[name][gene[rows]] + PSEUDOCOUNT)
                               / (cpm["NTC"][gene[rows]] + PSEUDOCOUNT))
    res["log2_fold_change"] = log2fc

    # Ensembl IDs that are not GENCODE v43 genes, and the ID the same gene carries there.
    updates = pd.read_csv(ID_UPDATES, sep="\t").set_index("ensembl_id").gencode_v43_id
    res["target_gene"] = res["target_gene"].replace(updates)

    res.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(res)} rows to {args.out}: {(res.p_val_adj < Q_THRESH).sum()} at q < {Q_THRESH}, "
          f"{np.isnan(log2fc).sum()} without log2_fold_change")


if __name__ == "__main__":
    main()
