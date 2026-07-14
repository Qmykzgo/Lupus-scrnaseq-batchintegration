"""
Smoke tests on small synthetic data. These do NOT touch the real GEO
download or run the full scanpy/scVI pipeline at scale — they exercise the
pure logic in each script (indexing, aggregation, stats) against a tiny
fake AnnData object so that import errors, wrong column names, and
off-by-one bugs surface immediately, without needing the real 1.26M-cell
file or a GPU.

Run:
    pytest tests/ -v

Requires the full environment (environment.yml / requirements.txt) to be
installed, since scanpy/anndata are used to build the fixture.
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import config  # noqa: E402
from utils import qc_metrics, resolve_gene_panel  # noqa: E402


# ---------------------------------------------------------------------------
# Fixture: a tiny synthetic version of the real object's structure
# ---------------------------------------------------------------------------
@pytest.fixture
def synthetic_adata():
    ad = pytest.importorskip("anndata")
    sc = pytest.importorskip("scanpy")
    rng = np.random.default_rng(0)

    n_cells = 400
    filler_genes = [f"GENE{i}" for i in range(40)]
    gene_names = filler_genes + config.CORE_ISG_PANEL  # ensure ISG panel is present
    n_genes = len(gene_names)

    counts = rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.float32)

    n_batches, n_individuals = 4, 20
    batches = rng.choice([f"pool{i}" for i in range(n_batches)], size=n_cells)
    cohorts = rng.choice(["cohort_A", "cohort_B"], size=n_cells)
    individuals = rng.choice([f"ind{i}" for i in range(n_individuals)], size=n_cells)
    # tie disease status to individual (real biology: it's a per-person label,
    # not independently random per cell)
    ind_disease = {f"ind{i}": rng.choice(["SLE", "Healthy"]) for i in range(n_individuals)}
    disease = np.array([ind_disease[i] for i in individuals])
    celltype = rng.choice(config.CELLTYPES_COARSE, size=n_cells)

    obs = pd.DataFrame(
        {
            config.BATCH_KEY: batches,
            config.PROCESSING_COHORT_KEY: cohorts,
            config.INDIVIDUAL_KEY: individuals,
            config.SAMPLE_KEY: [f"{i}_{b}" for i, b in zip(individuals, batches)],
            config.DISEASE_KEY: disease,
            config.CELLTYPE_COARSE_KEY: celltype,
        },
        index=[f"cell{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=gene_names)

    adata = ad.AnnData(X=counts.copy(), obs=obs, var=var)
    adata.raw = adata.copy()  # raw counts, matching the real object's structure
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)  # .X now log-normalized, .raw.X still raw counts
    return adata


# ---------------------------------------------------------------------------
# utils.py
# ---------------------------------------------------------------------------
def test_resolve_gene_panel_finds_all_isgs(synthetic_adata):
    present = resolve_gene_panel(synthetic_adata, config.CORE_ISG_PANEL)
    assert set(present) == set(config.CORE_ISG_PANEL)


def test_resolve_gene_panel_raises_when_mostly_absent(synthetic_adata):
    with pytest.raises(ValueError):
        resolve_gene_panel(synthetic_adata, ["NOT_A_REAL_GENE_1", "NOT_A_REAL_GENE_2"])


def test_qc_metrics_adds_expected_columns(synthetic_adata):
    out = qc_metrics(synthetic_adata)
    assert "n_genes_by_counts" in out.obs.columns
    assert "total_counts" in out.obs.columns
    assert (out.obs["total_counts"] >= 0).all()


# ---------------------------------------------------------------------------
# 01_subsample_dev_set.py — stratified_indices
# ---------------------------------------------------------------------------
def test_stratified_indices_preserves_all_batches(synthetic_adata):
    mod = importlib.import_module("01_subsample_dev_set")
    idx = mod.stratified_indices(synthetic_adata.obs, n_target=100, seed=0)
    assert len(idx) > 0
    subset_obs = synthetic_adata.obs.iloc[idx]
    # every batch present in the input should still have at least one cell
    # in a reasonably-sized subsample (not a strict guarantee at extreme
    # subsampling ratios, but true here given the fixture's batch sizes)
    assert set(subset_obs[config.BATCH_KEY]) == set(synthetic_adata.obs[config.BATCH_KEY])


def test_stratified_indices_respects_target_size_order_of_magnitude(synthetic_adata):
    mod = importlib.import_module("01_subsample_dev_set")
    idx = mod.stratified_indices(synthetic_adata.obs, n_target=100, seed=0)
    # stratified sampling won't hit the target exactly, but shouldn't wildly overshoot
    assert 20 <= len(idx) <= len(synthetic_adata)


# ---------------------------------------------------------------------------
# 02_qc_and_preprocessing.py — raw_counts_view
# ---------------------------------------------------------------------------
def test_raw_counts_view_recovers_counts_not_lognorm(synthetic_adata):
    mod = importlib.import_module("02_qc_and_preprocessing")
    raw_view = mod.raw_counts_view(synthetic_adata)
    # log-normalized values are (almost) never exact non-negative integers;
    # raw counts should be. This catches the classic bug of accidentally
    # using .X instead of .raw.X.
    X = raw_view.X.toarray() if hasattr(raw_view.X, "toarray") else raw_view.X
    assert np.allclose(X, np.round(X)), "raw_counts_view should return integer counts, not log-normalized values"
    assert raw_view.n_obs == synthetic_adata.n_obs


# ---------------------------------------------------------------------------
# 06_isg_downstream_analysis.py — cramers_v, build_pseudobulk, scoring
# ---------------------------------------------------------------------------
def test_cramers_v_zero_when_independent():
    mod = importlib.import_module("06_isg_downstream_analysis")
    # perfectly balanced contingency table -> independence -> V close to 0
    balanced = pd.DataFrame({"SLE": [50, 50], "Healthy": [50, 50]}, index=["cohort_A", "cohort_B"])
    v = mod.cramers_v(balanced)
    assert v < 0.05


def test_cramers_v_high_when_fully_confounded():
    mod = importlib.import_module("06_isg_downstream_analysis")
    confounded = pd.DataFrame({"SLE": [100, 0], "Healthy": [0, 100]}, index=["cohort_A", "cohort_B"])
    v = mod.cramers_v(confounded)
    assert v > 0.9


def test_score_isg_signature_adds_column(synthetic_adata):
    mod = importlib.import_module("06_isg_downstream_analysis")
    out = mod.score_isg_signature(synthetic_adata, logger=mod.logger)
    assert "ISG_score" in out.obs.columns
    assert out.obs["ISG_score"].notna().all()


def test_build_pseudobulk_shapes(synthetic_adata):
    mod = importlib.import_module("06_isg_downstream_analysis")
    synthetic_adata.layers["counts"] = synthetic_adata.raw.X.copy()
    cm_mask = synthetic_adata.obs[config.CELLTYPE_COARSE_KEY] == "cM"
    if cm_mask.sum() == 0:
        pytest.skip("Synthetic fixture happened to draw zero cM cells this run; re-seed if this recurs.")
    adata_cm = synthetic_adata[cm_mask].copy()
    pb_counts, meta = mod.build_pseudobulk(adata_cm, logger=mod.logger)
    assert pb_counts.shape[0] == meta.shape[0]
    assert set(meta.index) == set(pb_counts.index)
    assert (pb_counts.values >= 0).all()


def test_cramers_v_handles_single_category():
    mod = importlib.import_module("06_isg_downstream_analysis")
    # Single row/column contingency table should return 0.0 (no crash)
    single_col = pd.DataFrame({"SLE": [50, 50]}, index=["cohort_A", "cohort_B"])
    v = mod.cramers_v(single_col)
    assert v == 0.0

    single_row = pd.DataFrame({"SLE": [50], "Healthy": [50]}, index=["cohort_A"])
    v = mod.cramers_v(single_row)
    assert v == 0.0


def test_stratified_indices_extreme_subsampling(synthetic_adata):
    mod = importlib.import_module("01_subsample_dev_set")
    # Even with an extremely low target size, every batch should still have at least 1 cell
    idx = mod.stratified_indices(synthetic_adata.obs, n_target=5, seed=0)
    assert len(idx) > 0
    subset_obs = synthetic_adata.obs.iloc[idx]
    assert set(subset_obs[config.BATCH_KEY]) == set(synthetic_adata.obs[config.BATCH_KEY])


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
