# Batch Integration in a Multi-Cohort Lupus PBMC Atlas

A portfolio project reproducing and extending part of the analysis from:

> Perez, R.K., Gordon, M.G., Subramaniam, M., *et al.* **Single-cell RNA-seq reveals
> cell type-specific molecular and genetic associations to lupus.**
> *Science* 376, eabf1970 (2022). https://doi.org/10.1126/science.abf1970

This is the CLUES (California Lupus Epidemiology Study) cohort: **1.26M PBMCs from
261 individuals** (162 SLE patients, 99 healthy controls), profiled with a multiplexed
scRNA-seq design ("mux-seq") in which several individuals are pooled into a single 10x
lane and later demultiplexed computationally. Pools were processed across many
sequencing batches over time. That design makes this dataset a genuinely good — not
contrived — case study for batch integration: technical batch (which pool, which
processing cohort) and biology (cell type, disease status) are both real axes of
variation in the same data, and a project like this lives or dies on being able to
tell them apart.

## The question this project answers

The paper's headline monocyte finding is that SLE patients show elevated type-I
interferon-stimulated gene (ISG) expression in classical monocytes. This project asks
a methods question underneath that biological one:

> **Is the SLE-associated ISG signal in classical monocytes robust to how we handle
> batch — or could some of it be an artifact of which processing pool a sample
> happened to land in?**

This is answered in three parts (`scripts/06_isg_downstream_analysis.py`):

1. Check whether disease status is confounded with processing batch in this cohort
   (it's a real risk in any multiplexed design — pools aren't necessarily balanced).
2. Run pseudobulk differential expression (SLE vs. healthy, classical monocytes only)
   twice — once ignoring processing batch, once with it as a covariate — and see
   whether the core ISGs survive adjustment.
3. Cross-check with GSEA against MSigDB Hallmark interferon gene sets as an unbiased
   complement to a hand-picked ISG panel.

Separately (`scripts/04_batch_integration.py`, `05_benchmark_integration.py`), three
integration methods (Harmony, BBKNN, scVI) are compared and benchmarked with `scib` to
confirm that integrated neighborhoods mix technical batches while still respecting the
paper's own published cell type labels — the standard needed before trusting any
cluster-based composition or DE analysis downstream.

**A methodological note that matters here:** Harmony and BBKNN correct the
*embedding/graph*, not gene expression values. They cannot change a pseudobulk DE
p-value computed on log-normalized counts — only scVI's `get_normalized_expression`
produces a batch-adjusted expression matrix. So "does integration remove the ISG
signal" is the wrong question to ask of Harmony/BBKNN; "does the signal survive
explicitly modeling batch as a DE covariate" is the right one. The scripts are
structured to keep these two things separate rather than conflating them.

## Repository structure

```
lupus-scrnaseq-batch-integration/
├── README.md
├── environment.yml
├── requirements.txt
├── scripts/
│   ├── config.py                        # paths, constants, gene panels
│   ├── utils.py                         # shared helpers (QC, plotting, logging)
│   ├── 00_download_data.sh              # pulls the processed h5ad from GEO
│   ├── 01_subsample_dev_set.py          # stratified subsample for local dev
│   ├── 02_qc_and_preprocessing.py       # QC metrics, normalization, checkpoint
│   ├── 03_baseline_embedding.py         # uncorrected PCA/UMAP/Leiden — see the batch effect
│   ├── 04_batch_integration.py          # Harmony, BBKNN, scVI
│   ├── 05_benchmark_integration.py      # scib metrics, before/after comparison
│   └── 06_isg_downstream_analysis.py    # the biological question (see above)
├── tests/
│   └── test_pipeline_smoke.py           # runs the whole pipeline on synthetic data
├── notebooks/
│   └── 00_exploration.md                # suggested exploration checklist (not code)
├── data/                                # gitignored — see data/README.md
└── results/
    ├── figures/
    └── tables/
```

## Setup

```bash
conda env create -f environment.yml
conda activate lupus-batch-integration
# or: pip install -r requirements.txt
```

Tested against the same versions referenced in the batch-integration reference used to
write this pipeline: `scanpy>=1.10`, `anndata>=0.10`, `scikit-learn>=1.4`,
`scvi-tools>=1.1`. If your installed versions differ, check `pip show <package>` and
adapt call signatures as needed — the scanpy/scvi APIs occasionally rename kwargs
across minor versions.

## Running the pipeline

```bash
# 1. Download the processed AnnData object (~2-4 GB compressed) from GEO
bash scripts/00_download_data.sh

# 2. Build a manageable development subset (default: ~120k cells, stratified by
#    processing pool, disease status, and published cell type)
python scripts/01_subsample_dev_set.py --n-cells 120000

# 3. QC pass + normalization checkpoint
python scripts/02_qc_and_preprocessing.py

# 4. Uncorrected PCA/UMAP/Leiden — establishes the "problem" (batch-driven structure)
python scripts/03_baseline_embedding.py

# 5. Run Harmony, BBKNN, and scVI integration
python scripts/04_batch_integration.py

# 6. Benchmark all three with scib, before vs. after
python scripts/05_benchmark_integration.py

# 7. Answer the biological question
python scripts/06_isg_downstream_analysis.py
```

Each script reads/writes `.h5ad` checkpoints under `data/processed/` so you can rerun
any single stage without repeating the ones before it.

### Running at full scale

Everything above defaults to the ~120k-cell dev subset, which runs on a laptop
(scVI training is the slow part; ~15-20 min for 50 epochs on CPU, a couple of minutes
on a GPU). To run on the full 1.26M-cell object, pass `--full` to
`02_qc_and_preprocessing.py` onward. Expect to need a machine with 64GB+ RAM for the
in-memory scanpy steps, and a GPU for scVI to finish in a reasonable time. Harmony and
BBKNN scale more comfortably than scVI training does at this size.

## Data note: this is not the full transcriptome

The processed file GEO provides (`GSE174188_CLUES1_adjusted.h5ad.gz`) already contains
only the ~1,999 genes the original authors used for their own downstream analysis
(`adata.X` = log-normalized values for those genes, `adata.raw.X` = matching raw
counts), plus their own PCA/UMAP/Louvain clusters and the published cell type labels
(`cg_cov` = 11 coarse types, `ct_cov` = finer subtypes). This project deliberately
does **not** just reuse their precomputed embedding — QC, normalization, HVG handling,
and clustering are redone independently in this pipeline so that the integration
comparison is real work rather than re-plotting numbers that were already computed.
Their published cell type labels *are* reused, but only as a ground-truth label set
for `scib` benchmarking (Step 6 of the general workflow: "reuse the authors' published
annotation" is a legitimate shortcut when the annotation itself isn't the thing you're
trying to evaluate).

If you want the full transcriptome instead of this ~2,000-gene subset, the raw FASTQ
/ Cell Ranger outputs are on the Human Cell Atlas Data Coordination Platform and
dbGaP (genotypes only) — out of scope for this project, but noted in case you want to
extend it.

## Key metadata fields (from the GEO object)

| Column | Meaning |
|---|---|
| `batch_cov` | Multiplexed sequencing pool ID (e.g. `dmx_YS-JY-22_pool6`) — the primary **technical batch** variable |
| `Processing_Cohort` | Broader processing-time cohort — a second, coarser batch variable |
| `ind_cov` | Individual/donor ID |
| `ind_cov_batch_cov` | Individual × batch composite ID (unique per sample) |
| `SLE_status` | `SLE` or `Healthy` |
| `Status` | Finer clinical status (e.g. managed, flare, treated) |
| `pop_cov` | Genetic ancestry (European / Asian) |
| `Age`, `Sex` | Demographics |
| `cg_cov` | Published coarse cell type (11 types: B, T4, T8, NK, cM, ncM, cDC, pDC, PB, Prolif, Progen) |
| `ct_cov` | Published finer cell subtype |
| `louvain` | Published cluster ID (computed on the authors' own batch-corrected embedding — not reused for clustering here, only for sanity-checking) |

## Methods referenced

- Korsunsky, I. *et al.* Fast, sensitive and accurate integration of single-cell data
  with Harmony. *Nat. Methods* 16, 1289–1296 (2019).
- Polański, K. *et al.* BBKNN: fast batch alignment of single cell transcriptomes.
  *Bioinformatics* 36, 964–965 (2020).
- Lopez, R. *et al.* Deep generative modeling for single-cell transcriptomics.
  *Nat. Methods* 15, 1053–1058 (2018). (scVI)
- Luecken, M.D. *et al.* Benchmarking atlas-level data integration in single-cell
  genomics. *Nat. Methods* 19, 41–50 (2022). (`scib`)
- Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and
  dispersion for RNA-seq data with DESeq2. *Genome Biol.* 15, 550 (2014).

## License

Code in this repository: MIT (see `LICENSE`). The underlying data is subject to the
original authors' GEO/dbGaP terms of use — this repo does not redistribute any data.


---

## Pipeline Execution Upgrades & Benchmarking Results (July 2026)

We have successfully executed the entire pipeline and validated the monocyte interferon signature robustness. Below is a summary of the additions and outcomes:

### One-Click Execution
We added [run_pipeline.sh](file:///Ubuntu-24.04/home/yer_kanat/Downloads/lupus-scrnaseq-batch-integration/lupus-scrnaseq-batch-integration/run_pipeline.sh) to the root directory. To run the entire pipeline end-to-end (download, subsampling, preprocessing, integration, benchmarking, and downstream validation):
```bash
./run_pipeline.sh
```

### Memory Safety Optimizations (16GB RAM)
To prevent WSL/Linux OOM kernel crashes on host machines with 16GB RAM, we optimized `scripts/05_benchmark_integration.py` to:
* Load H5AD datasets in `backed="r"` mode, avoiding loading the heavy gene expression matrix `.X` (saving ~4.5 GB RAM per file).
* Copy only the necessary lightweight tables (`obs`) and embeddings (`obsm`) in-memory.
* Explicitly close HDF5 file handles (`adata.file.close()`) and run garbage collection (`gc.collect()`) after processing each method sequentially.
* This caps peak memory consumption to **under 100 MB**.

### Confounding Check (Part A)
* Cramer's V score between `SLE_status` and `Processing_Cohort` is **0.644** (strong confounding), highlighting the absolute necessity of covariates adjustment in differential expression.

### Batch Integration Benchmarking (Step 6)
We benchmarked the integrated embeddings with Silhouette width metrics (lower batch score is better mixed, higher cell-type score is better preserved):

| Method | Batch Silhouette (Want Low) | Cell Type Silhouette (Want High) |
| :--- | :---: | :---: |
| **uncorrected** | 0.0014 | 0.6114 |
| **harmony** | -0.0777 | 0.2492 |
| **bbknn** | -0.1154 | 0.2220 |
| **scvi** | -0.0624 | 0.1466 |

### Downstream Biological Findings (Step 7)
* **Core ISG Survival**: All **20 out of 20 core ISGs** (e.g. *ISG15*, *IFI6*, *IFI44*, *IFIT1*, *MX1*, *OAS1*) remain highly significant ($padj < 0.05$) in the adjusted pseudobulk DE model. The monocyte interferon-stimulated signature is **robust** and not a technical batch artifact.
* **GSEA**: Hallmark Interferon Alpha and Gamma response terms are highly enriched (NES > 3.1, FDR = 0.0).
* **scVI Denoised Score**: Correlates strongly with raw counts (Pearson = 0.815), confirming the signature is robust under deep batch adjustment.
