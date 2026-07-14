# data/

This directory is gitignored — nothing here is committed. Run the pipeline to
populate it.

```
data/
├── raw/
│   └── GSE174188_CLUES1_adjusted.h5ad     # from scripts/00_download_data.sh
└── processed/
    ├── lupus_dev_subset.h5ad              # from 01_subsample_dev_set.py
    ├── lupus_qc.h5ad                      # from 02_qc_and_preprocessing.py
    ├── lupus_baseline_embedding.h5ad      # from 03_baseline_embedding.py
    ├── lupus_harmony.h5ad                 # from 04_batch_integration.py
    ├── lupus_bbknn.h5ad
    └── lupus_scvi.h5ad
```

`raw/` should never be edited in place; every downstream script reads from
`processed/` and writes a new checkpoint back to `processed/`, so you can always
delete everything under `processed/` and rebuild it from `raw/` without
re-downloading.
