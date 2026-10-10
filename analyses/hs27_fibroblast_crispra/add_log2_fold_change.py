import argparse

import anndata as ad
import numpy as np
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser()
    # Output of `generate_element_level_results.py --padj_thresh 2`: every tested pair.
    p.add_argument('--results', default='element_level_results_unfiltered.tsv')
    # Count-scale regression from TFs_CRISPRa "Step 4 - Simplified regressions ... Hs27 experiment":
    # X holds OLS coefficients in UMIs per cell (guides x genes), with the baseline in row 'control'.
    p.add_argument('--count_coefs',
                   default='20240310_fibroblast_final_low_UMI_count_regressions_coef_mean_pop.h5ad')
    p.add_argument('--pseudocount_cpm', type=float, default=1)
    p.add_argument('--out', default='results/element_level_results_all_pairs.tsv')
    return p.parse_args()

def main():
    args = parse_args()

    # Read as text so every existing column is written back unchanged.
    res = pd.read_csv(args.results, sep='\t', dtype=str)

    adata = ad.read_h5ad(args.count_coefs, backed='r')
    guides = res.guide_id.unique().tolist()
    genes = res.target_gene.unique().tolist()
    # Control UMIs per cell summed over every gene, to express the pseudocount in CPM.
    control_total = np.clip(adata['control'].to_memory().X, 0, None).sum()
    sub = adata[guides + ['control'], genes].to_memory()
    coefs = pd.DataFrame(sub.X, index=sub.obs_names, columns=sub.var_names, dtype=np.float64)

    row = coefs.index.get_indexer(res.guide_id)
    col = coefs.columns.get_indexer(res.target_gene)
    control = coefs.loc['control'].values[col]
    # Modelled expression with the guide = baseline + guide coefficient.
    perturbed = control + coefs.values[row, col]

    # The additive model can predict negative expression for repressed genes with a low baseline;
    # expression cannot be negative, so clip at zero, then stabilise the ratio with a pseudocount.
    assert (control > 0).all()
    clipped = perturbed < 0
    pseudocount = args.pseudocount_cpm * control_total / 1e6
    res['log2_fold_change'] = np.log2((np.maximum(perturbed, 0) + pseudocount)
                                      / (control + pseudocount))

    res.to_csv(args.out, sep='\t', index=False)
    print(f"Wrote {len(res)} rows to {args.out}: pseudocount {args.pseudocount_cpm} CPM = "
          f"{pseudocount:.4f} UMIs per cell, {clipped.sum()} rows clipped at zero expression")

if __name__ == '__main__':
    main()
