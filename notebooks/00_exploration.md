# Exploration checklist

This project keeps its actual logic in `scripts/*.py` (importable, testable,
runnable end-to-end) rather than in notebooks, so that `tests/` can exercise
it directly. If you want an interactive notebook for exploring intermediate
checkpoints — which is genuinely useful for eyeballing plots while iterating
on thresholds — here's a suggested skeleton. Paste this into a fresh
`00_exploration.ipynb` after running at least through
`02_qc_and_preprocessing.py`:

```python
import sys
sys.path.insert(0, "../scripts")
import config
import scanpy as sc

adata = sc.read_h5ad(config.CKPT_QC)
adata

# Sanity-check the metadata fields the whole pipeline depends on
adata.obs[[config.BATCH_KEY, config.PROCESSING_COHORT_KEY,
           config.DISEASE_KEY, config.CELLTYPE_COARSE_KEY]].describe(include="all")

# How many cells per processing pool? (motivates why script 01 keeps every
# batch rather than uniform-random subsampling)
adata.obs[config.BATCH_KEY].value_counts().plot(kind="bar", figsize=(14, 4))

# QC distributions, per the workflow guide's convention of looking before
# choosing thresholds rather than copy-pasting them
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts"],
             groupby=config.PROCESSING_COHORT_KEY, rotation=45)
```

Once you've run `03_baseline_embedding.py` and `04_batch_integration.py`,
the most useful thing to look at interactively is the side-by-side UMAP grid
that `05_benchmark_integration.py` already saves to
`results/figures/05_integration_comparison_umaps.png` — start there before
re-deriving your own version in a notebook.

## Things worth eyeballing manually (not fully automated by the scripts)

- Do any single mux-seq pools look like outliers in the uncorrected UMAP
  (script 03) — e.g. one pool forming its own isolated island? That's a
  candidate for a closer look at that pool's QC metrics specifically.
- After integration (script 04/05): pick 2-3 pools and confirm their cells
  are now interleaved with other pools' cells of the same `cg_cov` type,
  not just globally "mixed" in aggregate.
- In the ISG analysis (script 06): plot `ISG_score` per individual, colored
  by `Processing_Cohort`, sorted by `SLE_status` — a quick visual gut-check
  before trusting the pseudobulk DE numbers.
