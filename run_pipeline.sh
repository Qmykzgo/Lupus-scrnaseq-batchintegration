#!/bin/bash
# Stop on any error
set -e

echo "========================================================================="
echo "       LUPUS SCRNA-SEQ BATCH INTEGRATION & ISG VALIDATION PIPELINE"
echo "========================================================================="
echo ""

# 1. Download data
echo ">>> [1/7] Downloading GSE174188 raw dataset from CELLxGENE..."
bash scripts/00_download_data.sh
echo ""

# 2. Subsample dev subset (30k cells for memory safety on 16GB RAM hosts)
echo ">>> [2/7] Subsampling development subset (30,000 cells)..."
python scripts/01_subsample_dev_set.py --n-cells 30000
echo ""

# 3. QC & Normalization
echo ">>> [3/7] Running Quality Control and normalization..."
python scripts/02_qc_and_preprocessing.py
echo ""

# 4. Baseline Embedding
echo ">>> [4/7] Generating uncorrected baseline embedding..."
python scripts/03_baseline_embedding.py
echo ""

# 5. Batch Integration (Harmony, BBKNN, scVI with GPU acceleration)
echo ">>> [5/7] Running batch integration methods (Harmony, BBKNN, scVI)..."
python scripts/04_batch_integration.py --scvi-gpu --scvi-epochs 50
echo ""

# 6. Integration Benchmarking
echo ">>> [6/7] Running integration benchmarking (silhouette metrics)..."
python scripts/05_benchmark_integration.py
echo ""

# 7. Downstream ISG Analysis
echo ">>> [7/7] Running downstream biological validation & scVI denoising..."
python scripts/06_isg_downstream_analysis.py --with-scvi-denoised
echo ""

echo "========================================================================="
echo "                      PIPELINE COMPLETED SUCCESSFULLY"
echo "========================================================================="
echo "Results and figures are saved under results/figures/ and results/tables/"
