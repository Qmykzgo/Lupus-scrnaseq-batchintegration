"""
Run three integration methods on `batch_cov` (the mux-seq pool ID) and
compute a corrected embedding + UMAP + Leiden clustering for each:

  - Harmony  (scanpy.external.pp.harmony_integrate) — fast, linear, corrects
    the PCA embedding.
  - BBKNN    (bbknn) — graph-based, corrects the kNN graph directly.
  - scVI     (scvi-tools) — deep generative model, learns a batch-invariant
    latent space from raw counts.

Each method writes its own checkpoint so script 05 can benchmark all three
against the shared uncorrected baseline from script 03. No single method is
assumed "correct" going in — see README and batch_integration reference for
why that comparison, rather than a single default choice, is the point.

Usage:
    python scripts/04_batch_integration.py                       # all three
    python scripts/04_batch_integration.py --methods harmony bbknn
    python scripts/04_batch_integration.py --scvi-epochs 100 --scvi-gpu
"""
import argparse

import config
import scanpy as sc
from utils import get_logger, require_file, save_checkpoint

logger = get_logger("04_integration")


def run_harmony(adata):
    import harmonypy as hm

    logger.info("Running Harmony (key=%s)...", config.BATCH_KEY)
    adata_h = adata.copy()
    
    pcs = adata_h.obsm["X_pca"]
    ho = hm.run_harmony(pcs, adata_h.obs, config.BATCH_KEY, random_state=config.RANDOM_STATE)
    
    # Z_corr is traditionally (n_pcs, n_cells), so Z_corr.T is (n_cells, n_pcs)
    corrected = ho.Z_corr.T
    if corrected.shape != (adata_h.n_obs, pcs.shape[1]):
        corrected = ho.Z_corr
        
    adata_h.obsm["X_pca_harmony"] = corrected
    
    sc.pp.neighbors(adata_h, use_rep="X_pca_harmony", n_neighbors=config.N_NEIGHBORS,
                     random_state=config.RANDOM_STATE)
    sc.tl.umap(adata_h, random_state=config.RANDOM_STATE)
    sc.tl.leiden(adata_h, resolution=config.LEIDEN_RESOLUTION,
                 random_state=config.RANDOM_STATE, key_added="leiden_harmony")
    return adata_h


def run_bbknn(adata):
    import bbknn

    logger.info("Running BBKNN (batch_key=%s)...", config.BATCH_KEY)
    adata_b = adata.copy()
    bbknn.bbknn(adata_b, batch_key=config.BATCH_KEY, n_pcs=config.N_PCS)
    sc.tl.umap(adata_b, random_state=config.RANDOM_STATE)
    sc.tl.leiden(adata_b, resolution=config.LEIDEN_RESOLUTION,
                 random_state=config.RANDOM_STATE, key_added="leiden_bbknn")
    return adata_b


def run_scvi(adata, max_epochs: int, use_gpu: bool):
    import scvi

    logger.info("Running scVI (batch_key=%s, max_epochs=%d)...", config.BATCH_KEY, max_epochs)
    if "counts" not in adata.layers:
        raise KeyError("adata.layers['counts'] not found — scVI needs raw counts. "
                        "Rerun scripts/02_qc_and_preprocessing.py.")
    adata_s = adata.copy()
    scvi.model.SCVI.setup_anndata(adata_s, batch_key=config.BATCH_KEY, layer="counts")
    model = scvi.model.SCVI(adata_s, n_latent=30, n_layers=2)
    accelerator = "gpu" if use_gpu else "cpu"
    model.train(max_epochs=max_epochs, early_stopping=True, accelerator=accelerator)
    
    # Save the trained model
    model_dir = config.DATA_PROCESSED / "scvi_model"
    model.save(model_dir, overwrite=True)
    logger.info("Saved scVI model to %s", model_dir)

    adata_s.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata_s, use_rep="X_scVI", n_neighbors=config.N_NEIGHBORS,
                     random_state=config.RANDOM_STATE)
    sc.tl.umap(adata_s, random_state=config.RANDOM_STATE)
    sc.tl.leiden(adata_s, resolution=config.LEIDEN_RESOLUTION,
                 random_state=config.RANDOM_STATE, key_added="leiden_scvi")
    return adata_s


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methods", nargs="+", default=["harmony", "bbknn", "scvi"],
                         choices=["harmony", "bbknn", "scvi"])
    parser.add_argument("--scvi-epochs", type=int, default=50,
                         help="~15-20 min for 50 epochs on CPU at dev-subset scale; "
                              "increase for production runs, especially with a GPU.")
    parser.add_argument("--scvi-gpu", action="store_true")
    args = parser.parse_args()

    require_file(config.CKPT_BASELINE, hint="Run scripts/03_baseline_embedding.py first.")
    adata = sc.read_h5ad(config.CKPT_BASELINE)
    logger.info("Loaded baseline checkpoint: %d cells x %d genes", adata.n_obs, adata.n_vars)

    if "harmony" in args.methods:
        save_checkpoint(run_harmony(adata), config.CKPT_HARMONY, logger)

    if "bbknn" in args.methods:
        save_checkpoint(run_bbknn(adata), config.CKPT_BBKNN, logger)

    if "scvi" in args.methods:
        save_checkpoint(run_scvi(adata, args.scvi_epochs, args.scvi_gpu),
                         config.CKPT_SCVI, logger)

    logger.info("Done. Run scripts/05_benchmark_integration.py to compare methods.")


if __name__ == "__main__":
    main()
