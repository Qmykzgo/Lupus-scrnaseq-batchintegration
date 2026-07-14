"""
QC + normalization, redone from raw counts rather than trusting the
pre-normalized `adata.X` shipped in the GEO object.

The downloaded object stores log-normalized values for ~1,999 genes in `.X`
and the matching raw counts in `.raw.X` (see README "Data note"). This script
starts from `.raw.X` so that normalization, HVG selection, and QC filtering
are this project's own work rather than a re-display of the original
authors' numbers — the point of a batch-integration project is to make
choices about that pipeline, including ones that affect what gets fed into
integration downstream.

Usage:
    python scripts/02_qc_and_preprocessing.py
    python scripts/02_qc_and_preprocessing.py --full   # start from the raw GEO file
"""
import argparse

import config
import scanpy as sc
from utils import get_logger, qc_metrics, require_file, save_checkpoint

logger = get_logger("02_qc")


def load_input(full: bool):
    if full:
        require_file(config.RAW_H5AD, hint="Run scripts/00_download_data.sh first.")
        path = config.RAW_H5AD
    else:
        require_file(config.CKPT_DEV_SUBSET,
                      hint="Run scripts/01_subsample_dev_set.py first.")
        path = config.CKPT_DEV_SUBSET
    logger.info("Loading %s", path)
    return sc.read_h5ad(path)


def raw_counts_view(adata):
    """Rebuild an AnnData whose .X is raw counts, using .raw (which the GEO
    object populates) rather than exponentiating the log-normalized .X back
    out — undoing a log1p+normalize_total round trip loses information that
    plain raw counts don't have."""
    if adata.raw is None:
        logger.warning(".raw not present on this object — falling back to .X "
                        "as-is. Downstream normalization will be a no-op if "
                        "this is already normalized data.")
        return adata.copy()
    counts_adata = adata.raw.to_adata()
    counts_adata.obs = adata.obs.copy()
    # obsm/uns from the original (e.g. any published PCA/UMAP) aren't needed
    # here since we're deliberately recomputing our own downstream.
    return counts_adata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true",
                         help="Read directly from the full raw GEO file "
                              "instead of the dev subset checkpoint.")
    args = parser.parse_args()

    adata = load_input(args.full)
    logger.info("Input: %d cells x %d genes", adata.n_obs, adata.n_vars)

    adata = raw_counts_view(adata)

    # --- QC metrics -----------------------------------------------------
    adata = qc_metrics(adata, logger)
    logger.info(
        "QC summary before filtering:\n%s",
        adata.obs[["n_genes_by_counts", "total_counts"]].describe().to_string(),
    )

    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=config.MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=config.MIN_CELLS_PER_GENE)
    if "pct_counts_mt" in adata.obs.columns and adata.obs["pct_counts_mt"].max() > 0:
        adata = adata[adata.obs["pct_counts_mt"] < config.MAX_PCT_MT].copy()
    logger.info("Cells retained after filtering: %d / %d (%.1f%%)",
                adata.n_obs, n_before, 100 * adata.n_obs / n_before)

    # --- Preserve raw counts, then normalize -----------------------------
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # --- Batch-aware HVG selection ---------------------------------------
    # Flags a subset of these ~1,999 genes as highly_variable *within* each
    # batch before combining, so PCA isn't driven by genes that only look
    # variable because of technical batch structure. Full gene matrix is
    # kept (not subsetted) — script 06 needs the ISG panel genes regardless
    # of whether they end up flagged highly_variable.
    sc.pp.highly_variable_genes(adata, batch_key=config.BATCH_KEY)
    logger.info("Highly variable genes flagged: %d / %d",
                int(adata.var["highly_variable"].sum()), adata.n_vars)

    save_checkpoint(adata, config.CKPT_QC, logger)


if __name__ == "__main__":
    main()
