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

    loaded = {}
    for name, (path, embed_key) in METHOD_CHECKPOINTS.items():
        if not path.exists():
            logger.warning("Checkpoint for '%s' not found at %s — skipping. "
                            "Run scripts/04_batch_integration.py.", name, path)
            continue
        loaded[name] = (sc.read_h5ad(path), embed_key)
        logger.info("Loaded '%s': %d cells", name, loaded[name][0].n_obs)

    if "uncorrected" not in loaded:
        raise SystemExit("Need at least the uncorrected baseline checkpoint to compare against.")

    adata_uncorrected = loaded["uncorrected"][0]

    rows = []
    for name, (adata, embed_key) in loaded.items():
        if name == "uncorrected":
            continue
        logger.info("Benchmarking '%s' (embedding: %s)...", name, embed_key)

        scib_metrics = try_scib_metrics(adata_uncorrected, adata, embed_key)
        manual = manual_silhouette(adata, embed_key)

        row = {"method": name, **manual}
        if scib_metrics is not None:
            # scib returns a DataFrame/Series depending on version; normalize to dict
            row["scib_metrics"] = (
                scib_metrics.to_dict() if hasattr(scib_metrics, "to_dict") else dict(scib_metrics)
            )
        rows.append(row)

    summary = pd.DataFrame(rows).set_index("method")
    logger.info("Integration benchmark summary:\n%s", summary.to_string())
    save_table(summary, "05_integration_benchmark_summary", logger)

    # --- Side-by-side UMAP comparison -------------------------------------
    # Per the workflow guide: good mixing AND clean cell-type separation
    # together is the signature of good integration; mixing alone can mean
    # overcorrection, so both panels matter.
    import matplotlib.pyplot as plt

    methods_present = [m for m in ["uncorrected", "harmony", "bbknn", "scvi"] if m in loaded]
    fig, axes = plt.subplots(2, len(methods_present), figsize=(5 * len(methods_present), 9))
    for col, name in enumerate(methods_present):
        adata = loaded[name][0]
        sc.pl.umap(adata, color=config.BATCH_KEY, ax=axes[0, col], show=False,
                   title=f"{name}: batch_cov", legend_loc=None, size=3)
        sc.pl.umap(adata, color=config.CELLTYPE_COARSE_KEY, ax=axes[1, col], show=False,
                   title=f"{name}: cg_cov", legend_loc=None if col < len(methods_present) - 1 else "right margin",
                   size=3)
    fig.tight_layout()
    save_figure(fig, "05_integration_comparison_umaps", logger)


if __name__ == "__main__":
    main()
