"""
Uncorrected PCA -> kNN -> UMAP -> Leiden, with no batch correction applied.

This establishes the problem that scripts 04-05 then address: how much of
the clustering structure in this cohort is driven by which mux-seq pool a
cell came from, rather than by cell type or disease status. Produces the
"before" comparison point for the integration benchmarking in script 05.

Usage:
    python scripts/03_baseline_embedding.py
"""
import config
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score
from utils import get_logger, require_file, save_checkpoint, save_figure, save_table, setup_plot_style

logger = get_logger("03_baseline")


def main():
    require_file(config.CKPT_QC, hint="Run scripts/02_qc_and_preprocessing.py first.")
    adata = sc.read_h5ad(config.CKPT_QC)
    logger.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)

    sc.pp.scale(adata, max_value=10)  # PCA benefits from scaling; HVG flag from script 02 is respected
    sc.pp.pca(adata, n_comps=config.N_PCS, random_state=config.RANDOM_STATE)
    sc.pp.neighbors(adata, n_neighbors=config.N_NEIGHBORS, n_pcs=config.N_PCS,
                     random_state=config.RANDOM_STATE)
    sc.tl.umap(adata, random_state=config.RANDOM_STATE)
    sc.tl.leiden(adata, resolution=config.LEIDEN_RESOLUTION,
                 random_state=config.RANDOM_STATE, key_added="leiden_uncorrected")

    n_clusters = adata.obs["leiden_uncorrected"].nunique()
    logger.info("Uncorrected Leiden clustering: %d clusters", n_clusters)

    # --- Quantify batch- vs biology-driven structure ---------------------
    # Adjusted Rand Index between the uncorrected clustering and (a) technical
    # batch, (b) published cell type. A clustering that agrees more with
    # batch_cov than with cg_cov is a red flag that batch effects, not
    # biology, are dominating — exactly the failure mode integration should
    # fix. This is a quick diagnostic, not a substitute for the scib
    # benchmarking in script 05.
    ari_batch = adjusted_rand_score(adata.obs[config.BATCH_KEY], adata.obs["leiden_uncorrected"])
    ari_celltype = adjusted_rand_score(adata.obs[config.CELLTYPE_COARSE_KEY], adata.obs["leiden_uncorrected"])
    logger.info("ARI(uncorrected clusters, batch_cov)   = %.4f", ari_batch)
    logger.info("ARI(uncorrected clusters, cg_cov)      = %.4f", ari_celltype)

    save_table(
        pd.DataFrame({
            "comparison": ["leiden_uncorrected vs batch_cov", "leiden_uncorrected vs cg_cov"],
            "adjusted_rand_index": [ari_batch, ari_celltype],
        }),
        "baseline_ari_batch_vs_celltype", logger,
    )

    # --- Plots -------------------------------------------------------------
    setup_plot_style()
    fig = sc.pl.umap(
        adata,
        color=[config.BATCH_KEY, config.PROCESSING_COHORT_KEY,
               config.CELLTYPE_COARSE_KEY, config.DISEASE_KEY],
        ncols=2, wspace=0.4, show=False, return_fig=True,
        legend_fontsize=6,
    )
    save_figure(fig, "03_baseline_umap_uncorrected", logger)

    save_checkpoint(adata, config.CKPT_BASELINE, logger)


if __name__ == "__main__":
    main()
