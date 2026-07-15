"""
Benchmark the three integrated embeddings against the uncorrected baseline
using `scib`, with a manual sklearn-silhouette fallback in case the
installed `scib` version's API doesn't match (its function signatures have
moved around across releases — see the batch-integration reference's version
note).

Two things are measured for every method, because optimizing only one is a
known failure mode:
  - Batch mixing:      do cells from different pools now sit near each
                        other? (higher = better mixed)
  - Bio conservation:  do cells of different published cell types (cg_cov)
                        still separate? (higher = better preserved)
A method that maximizes mixing while destroying cell-type separation has
overcorrected — it hasn't "integrated," it's erased the biology.

Usage:
    python scripts/05_benchmark_integration.py
"""
import config
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_score
from utils import get_logger, save_figure, save_table, setup_plot_style

logger = get_logger("05_benchmark")

METHOD_CHECKPOINTS = {
    "uncorrected": (config.CKPT_BASELINE, "X_pca"),
    "harmony": (config.CKPT_HARMONY, "X_pca_harmony"),
    "bbknn": (config.CKPT_BBKNN, "X_pca"),   # BBKNN corrects the graph, not the embedding
    "scvi": (config.CKPT_SCVI, "X_scVI"),
}


def try_scib_metrics(adata_unintegrated, adata_integrated, embed_key: str):
    """Attempt the full scib metric suite; return None on any failure so the
    sklearn fallback below still gives a usable (if smaller) comparison."""
    try:
        import scib
        metrics = scib.metrics.metrics_fast(
            adata_unintegrated, adata_integrated,
            batch_key=config.BATCH_KEY, label_key=config.CELLTYPE_COARSE_KEY,
            embed=embed_key,
        )
        return metrics
    except Exception as exc:  # noqa: BLE001 - deliberately broad; this is a
        # version-compatibility fallback, not a correctness bug to hide.
        logger.warning("scib.metrics.metrics_fast failed (%s: %s) — falling back "
                        "to a manual silhouette comparison. If you hit this, run "
                        "`pip show scib` and check scib.metrics' current function "
                        "signatures; the API has changed across releases.",
                        type(exc).__name__, exc)
        return None


def manual_silhouette(adata, embed_key: str):
    """Matches the pattern in the batch-integration reference: batch
    silhouette should be LOW (well mixed), cell-type silhouette should be
    HIGH (well separated). Subsamples for tractable runtime on large N."""
    n = adata.n_obs
    sample = adata if n <= 20_000 else sc.pp.subsample(adata, n_obs=20_000, copy=True)
    X = sample.obsm[embed_key]
    batch_sil = silhouette_score(X, sample.obs[config.BATCH_KEY])
    celltype_sil = silhouette_score(X, sample.obs[config.CELLTYPE_COARSE_KEY])
    return {"batch_silhouette (want low)": batch_sil,
            "celltype_silhouette (want high)": celltype_sil}


def main():
    setup_plot_style()
    import gc

    # 1. Load uncorrected baseline first in backed mode to avoid loading heavy .X dense matrix
    uncorrected_path, uncorrected_embed = METHOD_CHECKPOINTS["uncorrected"]
    if not uncorrected_path.exists():
        raise SystemExit(f"Need uncorrected baseline checkpoint at {uncorrected_path} to compare against.")
    
    logger.info("Loading baseline checkpoint in backed mode...")
    adata_backed = sc.read_h5ad(uncorrected_path, backed="r")
    adata_uncorrected = sc.AnnData(
        obs=adata_backed.obs.copy(),
        obsm={k: adata_backed.obsm[k].copy() for k in adata_backed.obsm.keys()}
    )
    adata_backed.file.close()
    del adata_backed
    logger.info("Loaded 'uncorrected': %d cells", adata_uncorrected.n_obs)

    # We will keep light versions of each method for final plotting to save RAM
    plot_adata = {
        "uncorrected": sc.AnnData(
            obs=adata_uncorrected.obs.copy(),
            obsm={"X_umap": adata_uncorrected.obsm["X_umap"].copy()}
        )
    }

    rows = []
    for name, (path, embed_key) in METHOD_CHECKPOINTS.items():
        if name == "uncorrected":
            continue
        if not path.exists():
            logger.warning("Checkpoint for '%s' not found at %s — skipping.", name, path)
            continue
        
        logger.info("Loading checkpoint for '%s' in backed mode...", name)
        adata_backed = sc.read_h5ad(path, backed="r")
        adata = sc.AnnData(
            obs=adata_backed.obs.copy(),
            obsm={k: adata_backed.obsm[k].copy() for k in adata_backed.obsm.keys()}
        )
        adata_backed.file.close()
        del adata_backed
        logger.info("Loaded '%s': %d cells", name, adata.n_obs)

        logger.info("Benchmarking '%s' (embedding: %s)...", name, embed_key)
        scib_metrics = try_scib_metrics(adata_uncorrected, adata, embed_key)
        manual = manual_silhouette(adata, embed_key)

        row = {"method": name, **manual}
        if scib_metrics is not None:
            row["scib_metrics"] = (
                scib_metrics.to_dict() if hasattr(scib_metrics, "to_dict") else dict(scib_metrics)
            )
        rows.append(row)

        # Store light copy for plotting and free full memory
        plot_adata[name] = sc.AnnData(
            obs=adata.obs.copy(),
            obsm={"X_umap": adata.obsm["X_umap"].copy()}
        )
        
        del adata
        gc.collect()

    summary = pd.DataFrame(rows).set_index("method")
    logger.info("Integration benchmark summary:\n%s", summary.to_string())
    save_table(summary, "05_integration_benchmark_summary", logger)

    # --- Side-by-side UMAP comparison -------------------------------------
    import matplotlib.pyplot as plt

    methods_present = [m for m in ["uncorrected", "harmony", "bbknn", "scvi"] if m in plot_adata]
    fig, axes = plt.subplots(2, len(methods_present), figsize=(5 * len(methods_present), 9))
    for col, name in enumerate(methods_present):
        adata = plot_adata[name]
        sc.pl.umap(adata, color=config.BATCH_KEY, ax=axes[0, col], show=False,
                   title=f"{name}: batch_cov", legend_loc=None, size=3)
        sc.pl.umap(adata, color=config.CELLTYPE_COARSE_KEY, ax=axes[1, col], show=False,
                   title=f"{name}: cg_cov", legend_loc=None if col < len(methods_present) - 1 else "right margin",
                   size=3)
    fig.tight_layout()
    save_figure(fig, "05_integration_comparison_umaps", logger)


if __name__ == "__main__":
    main()
