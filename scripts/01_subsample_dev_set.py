"""
Build a manageable development subset from the full 1.26M-cell object.

Rather than a uniform random cell sample, this keeps EVERY processing batch
(`batch_cov`) but downsamples proportionally within each
(batch_cov x SLE_status x cg_cov) stratum. That preserves the thing this
project actually needs: every pool still represented, with realistic relative
proportions of disease status and cell type per pool, so downstream batch
integration is solving a smaller version of the real problem rather than an
easier synthetic one.

Uses AnnData's backed mode to select rows without loading the full
1.26M x 1,999 matrix into memory first — reading only what's needed.

Usage:
    python scripts/01_subsample_dev_set.py --n-cells 120000
    python scripts/01_subsample_dev_set.py --full     # skip subsampling entirely
"""
import argparse

import config
import numpy as np
import scanpy as sc
from utils import get_logger, require_file, save_checkpoint

logger = get_logger("01_subsample")


def stratified_indices(obs, n_target: int, seed: int) -> np.ndarray:
    strata_cols = [config.BATCH_KEY, config.DISEASE_KEY, config.CELLTYPE_COARSE_KEY]
    missing = [c for c in strata_cols if c not in obs.columns]
    if missing:
        raise KeyError(f"Expected stratification columns not found in .obs: {missing}")

    frac = min(1.0, n_target / len(obs))
    logger.info(
        "Stratified sampling at fraction=%.4f across %s (target ~%d of %d cells)",
        frac, strata_cols, n_target, len(obs),
    )
    
    def sample_group(g):
        if len(g) == 0:
            return g
        n = max(1, int(np.round(len(g) * frac)))
        n = min(len(g), n)
        return g.sample(n=n, random_state=seed)

    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning, message=".*DataFrameGroupBy.apply.*")
        sampled = obs.groupby(strata_cols, observed=True, group_keys=False).apply(sample_group)
        
    idx = obs.index.get_indexer(sampled.index)
    return np.sort(idx)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-cells", type=int, default=config.DEV_N_CELLS_DEFAULT,
                         help="Approximate target cell count for the dev subset.")
    parser.add_argument("--full", action="store_true",
                         help="Skip subsampling; just copy the full object forward "
                              "as the checkpoint (only recommended on a machine "
                              "with 64GB+ RAM and, ideally, a GPU for scVI).")
    parser.add_argument("--seed", type=int, default=config.DEV_SEED)
    args = parser.parse_args()

    require_file(
        config.RAW_H5AD,
        hint="Run scripts/00_download_data.sh first.",
    )

    logger.info("Opening %s in backed mode (not loading full matrix into memory)...",
                config.RAW_H5AD)
    adata = sc.read_h5ad(config.RAW_H5AD, backed="r")
    logger.info("Full object: %d cells x %d genes", adata.n_obs, adata.n_vars)

    # Map CELLxGENE columns to the author's original names expected by the pipeline
    obs = adata.obs.copy()
    mapped = False
    if 'library_uuid' in obs.columns and config.BATCH_KEY not in obs.columns:
        logger.info("Mapping CELLxGENE 'library_uuid' -> '%s'", config.BATCH_KEY)
        obs[config.BATCH_KEY] = obs['library_uuid']
        mapped = True
    if 'author_cell_type' in obs.columns and config.CELLTYPE_COARSE_KEY not in obs.columns:
        logger.info("Mapping CELLxGENE 'author_cell_type' -> '%s'", config.CELLTYPE_COARSE_KEY)
        obs[config.CELLTYPE_COARSE_KEY] = obs['author_cell_type']
        mapped = True
    if 'disease' in obs.columns and config.DISEASE_KEY not in obs.columns:
        logger.info("Mapping CELLxGENE 'disease' -> '%s'", config.DISEASE_KEY)
        obs[config.DISEASE_KEY] = obs['disease'].map({
            'systemic lupus erythematosus': 'SLE',
            'normal': 'Healthy'
        })
        mapped = True
    if config.SAMPLE_KEY not in obs.columns and config.INDIVIDUAL_KEY in obs.columns and config.BATCH_KEY in obs.columns:
        obs[config.SAMPLE_KEY] = obs[config.INDIVIDUAL_KEY].astype(str) + "_" + obs[config.BATCH_KEY].astype(str)
        mapped = True
        
    if mapped:
        adata.obs = obs

    for col in (config.BATCH_KEY, config.DISEASE_KEY, config.CELLTYPE_COARSE_KEY,
                config.INDIVIDUAL_KEY, config.PROCESSING_COHORT_KEY):
        if col not in adata.obs.columns:
            logger.warning("Expected column '%s' not found in .obs — check that "
                            "the downloaded file matches the version documented "
                            "in the README.", col)

    if args.full:
        logger.info("--full passed: loading entire object into memory (this is large).")
        subset = adata.to_memory()
    else:
        idx = stratified_indices(adata.obs, args.n_cells, args.seed)
        logger.info("Selected %d cells; reading just those rows from disk...", len(idx))
        subset = adata[idx].to_memory()

    logger.info("Dev subset: %d cells x %d genes, %d individuals, %d batches",
                subset.n_obs, subset.n_vars,
                subset.obs[config.INDIVIDUAL_KEY].nunique(),
                subset.obs[config.BATCH_KEY].nunique())
    logger.info("Disease status counts:\n%s",
                subset.obs[config.DISEASE_KEY].value_counts().to_string())

    save_checkpoint(subset, config.CKPT_DEV_SUBSET, logger)


if __name__ == "__main__":
    main()
