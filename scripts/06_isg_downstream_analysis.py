"""
The biological question this project is built around (see README):

    Is the SLE-associated interferon-stimulated gene (ISG) signal in
    classical monocytes robust to processing batch, or could some of it be
    a batch confound?

Four parts, run in order:
  A. Is SLE_status confounded with Processing_Cohort in this cohort at all?
     (quantified with Cramer's V on the contingency table)
  B. Per-cell ISG signature score in monocytes (cM/ncM), SLE vs. Healthy,
     visualized against both disease status and processing cohort.
  C. Pseudobulk differential expression in classical monocytes (cM) only,
     SLE vs. Healthy, run twice: naively (~ SLE_status) and adjusted
     (~ Processing_Cohort + SLE_status). Compares which core ISGs survive
     covariate adjustment.
  D. GSEA against MSigDB Hallmark interferon response gene sets on the
     adjusted model's ranked gene statistics, as an unbiased cross-check
     against the hand-picked ISG panel in config.py.

Important: this script deliberately works on `data/processed/lupus_qc.h5ad`
(expression values), not on any of the *_harmony/_bbknn/_scvi checkpoints —
Harmony and BBKNN correct the embedding/graph, not gene expression, so they
have nothing to contribute to a DE analysis on expression values. Only
scVI's `get_normalized_expression` produces a batch-adjusted expression
matrix; that comparison is offered as an optional extra at the end.

Usage:
    python scripts/06_isg_downstream_analysis.py
    python scripts/06_isg_downstream_analysis.py --with-scvi-denoised
"""
import argparse

import config
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import chi2_contingency
from utils import get_logger, require_file, resolve_gene_panel, save_figure, save_table, setup_plot_style

logger = get_logger("06_isg_analysis")


# ---------------------------------------------------------------------------
# Part A — confounding check
# ---------------------------------------------------------------------------
def cramers_v(contingency: pd.DataFrame) -> float:
    if contingency.shape[0] <= 1 or contingency.shape[1] <= 1:
        return 0.0
    chi2, _, _, _ = chi2_contingency(contingency)
    n = contingency.values.sum()
    if n == 0:
        return 0.0
    r, c = contingency.shape
    return float(np.sqrt((chi2 / n) / (min(r, c) - 1)))


def check_confounding(adata, logger):
    per_individual = (
        adata.obs[[config.INDIVIDUAL_KEY, config.DISEASE_KEY, config.PROCESSING_COHORT_KEY]]
        .drop_duplicates(subset=config.INDIVIDUAL_KEY)
    )
    ct = pd.crosstab(per_individual[config.PROCESSING_COHORT_KEY], per_individual[config.DISEASE_KEY])
    v = cramers_v(ct)
    logger.info("SLE_status x Processing_Cohort contingency table (by individual):\n%s", ct.to_string())
    logger.info("Cramer's V (association strength, 0=independent, 1=fully confounded) = %.3f", v)
    if v > 0.3:
        logger.warning("Moderate-to-strong association between disease status and processing "
                        "cohort in this cohort — the adjusted DE model in Part C matters, not "
                        "just as a formality.")
    save_table(ct, "06_confounding_contingency_table", logger)
    return v


# ---------------------------------------------------------------------------
# Part B — per-cell signature scoring
# ---------------------------------------------------------------------------
def score_isg_signature(adata, logger):
    panel = resolve_gene_panel(adata, config.CORE_ISG_PANEL, logger)
    sc.tl.score_genes(adata, gene_list=panel, score_name="ISG_score",
                       random_state=config.RANDOM_STATE)
    return adata


def plot_isg_by_group(adata_mono, logger):
    import matplotlib.pyplot as plt
    import seaborn as sns

    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    sns.boxplot(data=adata_mono.obs, x=config.DISEASE_KEY, y="ISG_score", ax=axes[0])
    axes[0].set_title("ISG score by disease status (cM + ncM)")
    sns.boxplot(data=adata_mono.obs, x=config.PROCESSING_COHORT_KEY, y="ISG_score", ax=axes[1])
    axes[1].set_title("ISG score by processing cohort (cM + ncM)")
    axes[1].tick_params(axis="x", rotation=45)
    fig.tight_layout()
    save_figure(fig, "06_isg_score_by_group", logger)


# ---------------------------------------------------------------------------
# Part C — pseudobulk DE, naive vs. batch-adjusted
# ---------------------------------------------------------------------------
def build_pseudobulk(adata_cm, logger):
    """Sum raw counts per individual for classical monocytes only. Uses
    `ind_cov` (individual), not `ind_cov_batch_cov`, as the replicate unit —
    the biological unit of replication is the person, per the workflow
    guide's pseudoreplication principle."""
    if "counts" not in adata_cm.layers:
        raise KeyError("layers['counts'] missing — rerun scripts/02_qc_and_preprocessing.py.")

    counts = pd.DataFrame(
        adata_cm.layers["counts"].toarray()
        if hasattr(adata_cm.layers["counts"], "toarray") else adata_cm.layers["counts"],
        index=adata_cm.obs_names, columns=adata_cm.var_names,
    )
    counts[config.INDIVIDUAL_KEY] = adata_cm.obs[config.INDIVIDUAL_KEY].values
    pb_counts = counts.groupby(config.INDIVIDUAL_KEY).sum()

    n_cells_per_ind = adata_cm.obs.groupby(config.INDIVIDUAL_KEY, observed=True).size()
    keep = n_cells_per_ind[n_cells_per_ind >= 10].index  # drop individuals with too few cM cells
    pb_counts = pb_counts.loc[pb_counts.index.intersection(keep)]

    meta = (
        adata_cm.obs[[config.INDIVIDUAL_KEY, config.DISEASE_KEY, config.PROCESSING_COHORT_KEY]]
        .drop_duplicates(subset=config.INDIVIDUAL_KEY)
        .set_index(config.INDIVIDUAL_KEY)
        .loc[pb_counts.index]
    )
    logger.info("Pseudobulk matrix: %d individuals (>=10 cM cells each) x %d genes",
                pb_counts.shape[0], pb_counts.shape[1])
    logger.info("Disease status counts in pseudobulk cohort:\n%s",
                meta[config.DISEASE_KEY].value_counts().to_string())
    return pb_counts, meta


def run_pseudobulk_de(counts_df, metadata_df, design_factors, contrast, logger, label):
    """pydeseq2's kwarg names have shifted across versions (design_factors
    vs. design=formula string). Try the modern API first and fall back to
    the older one rather than hard failing — matches the version-tolerance
    convention used elsewhere in this pipeline."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts_df = counts_df.round().astype(int)  # pydeseq2 requires integer counts
    try:
        dds = DeseqDataSet(counts=counts_df, metadata=metadata_df, design_factors=design_factors)
    except TypeError:
        logger.info("DeseqDataSet(design_factors=...) not accepted by this pydeseq2 version; "
                    "trying design=<formula string> instead. Check `help(DeseqDataSet)` if this "
                    "also fails.")
        formula = "~" + " + ".join(design_factors if isinstance(design_factors, list) else [design_factors])
        dds = DeseqDataSet(counts=counts_df, metadata=metadata_df, design=formula)

    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=contrast)
    stat_res.summary()
    results = stat_res.results_df.sort_values("padj")
    logger.info("[%s] Top 10 DE genes:\n%s", label, results.head(10).to_string())
    save_table(results, f"06_deseq2_{label}", logger)
    return results


def compare_isg_survival(results_naive, results_adjusted, isg_panel, logger):
    def sig_isgs(results):
        present = [g for g in isg_panel if g in results.index]
        sig = results.loc[present]
        return set(sig[sig["padj"] < 0.05].index)

    naive_sig = sig_isgs(results_naive)
    adjusted_sig = sig_isgs(results_adjusted)
    logger.info("Core ISGs significant (padj<0.05), naive model (~SLE_status): %s", sorted(naive_sig))
    logger.info("Core ISGs significant (padj<0.05), adjusted model (~Processing_Cohort + SLE_status): %s",
                sorted(adjusted_sig))
    lost = naive_sig - adjusted_sig
    if lost:
        logger.warning("ISGs significant naively but NOT after batch adjustment (possible batch "
                        "confound rather than genuine disease signal): %s", sorted(lost))
    else:
        logger.info("All naively-significant core ISGs survive batch adjustment — signal looks robust.")
    return pd.DataFrame({
        "gene": sorted(naive_sig | adjusted_sig),
    }).assign(
        significant_naive=lambda d: d["gene"].isin(naive_sig),
        significant_adjusted=lambda d: d["gene"].isin(adjusted_sig),
    )


# ---------------------------------------------------------------------------
# Part D — GSEA cross-check
# ---------------------------------------------------------------------------
def run_gsea(results_adjusted, logger):
    import gseapy as gp

    ranked = results_adjusted["stat"].dropna().sort_values(ascending=False)
    try:
        pre_res = gp.prerank(rnk=ranked, gene_sets=["MSigDB_Hallmark_2020"],
                              permutation_num=1000, outdir=None, seed=config.RANDOM_STATE)
    except Exception as exc:  # noqa: BLE001 — network access to Enrichr is required here;
        # fail informatively rather than crash the whole pipeline run.
        logger.warning("GSEA step failed (%s: %s). This step needs internet access to fetch "
                        "MSigDB_Hallmark_2020 from Enrichr at runtime — skipping.",
                        type(exc).__name__, exc)
        return None

    res2d = pre_res.res2d
    ifn_terms = res2d[res2d["Term"].str.contains("Interferon", case=False, na=False)]
    logger.info("Hallmark interferon-related GSEA terms:\n%s", ifn_terms.to_string())
    save_table(res2d, "06_gsea_hallmark_full", logger)
    save_table(ifn_terms, "06_gsea_hallmark_interferon_terms", logger)
    return ifn_terms


# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--with-scvi-denoised", action="store_true",
                         help="Also compare against scVI's model-based denoised expression "
                              "(the one integration method that DOES touch expression values).")
    args = parser.parse_args()

    require_file(config.CKPT_QC, hint="Run scripts/02_qc_and_preprocessing.py first.")
    adata = sc.read_h5ad(config.CKPT_QC)
    logger.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)

    logger.info("--- Part A: confounding check ---")
    check_confounding(adata, logger)

    logger.info("--- Part B: per-cell ISG signature scoring ---")
    adata = score_isg_signature(adata, logger)
    adata_mono = adata[adata.obs[config.CELLTYPE_COARSE_KEY].isin(config.MONOCYTE_TYPES)].copy()
    plot_isg_by_group(adata_mono, logger)

    logger.info("--- Part C: pseudobulk DE, classical monocytes only ---")
    adata_cm = adata[adata.obs[config.CELLTYPE_COARSE_KEY] == "cM"].copy()
    pb_counts, meta = build_pseudobulk(adata_cm, logger)

    results_naive = run_pseudobulk_de(
        pb_counts, meta, design_factors=[config.DISEASE_KEY],
        contrast=[config.DISEASE_KEY, "SLE", "Healthy"], logger=logger, label="naive",
    )
    results_adjusted = run_pseudobulk_de(
        pb_counts, meta, design_factors=[config.PROCESSING_COHORT_KEY, config.DISEASE_KEY],
        contrast=[config.DISEASE_KEY, "SLE", "Healthy"], logger=logger, label="adjusted",
    )
    isg_panel_present = resolve_gene_panel(adata_cm, config.CORE_ISG_PANEL, logger)
    survival_table = compare_isg_survival(results_naive, results_adjusted, isg_panel_present, logger)
    save_table(survival_table, "06_isg_survival_naive_vs_adjusted", logger)

    logger.info("--- Part D: GSEA cross-check against MSigDB Hallmark ---")
    run_gsea(results_adjusted, logger)

    if args.with_scvi_denoised:
        logger.info("--- Optional: scVI denoised-expression comparison ---")
        run_scvi_denoised_comparison(logger)

    logger.info("Done. See results/tables/06_* and results/figures/06_* for output.")


def run_scvi_denoised_comparison(logger):
    """scVI is the one integration method in this project that can produce a
    batch-adjusted expression matrix (via get_normalized_expression), unlike
    Harmony/BBKNN which only correct the embedding/graph. This loads the saved
    trained scVI model, extracts denoised normalized expression, and scores
    cells to see if the SLE monocyte signature persists."""
    import scvi
    import matplotlib.pyplot as plt
    import seaborn as sns

    model_dir = config.DATA_PROCESSED / "scvi_model"
    require_file(model_dir, hint="Run scripts/04_batch_integration.py --methods scvi first.")
    
    logger.info("Loading trained scVI model from %s...", model_dir)
    require_file(config.CKPT_QC, hint="Run scripts/02_qc_and_preprocessing.py first.")
    adata = sc.read_h5ad(config.CKPT_QC)
    
    # scvi needs the adata to match the setup
    model = scvi.model.SCVI.load(model_dir, adata=adata)
    logger.info("Extracting denoised normalized expression from scVI model...")
    denoised_expr = model.get_normalized_expression()
    
    # Create an AnnData with scVI denoised expression in .X
    adata_denoised = adata.copy()
    adata_denoised.X = denoised_expr
    
    logger.info("Scoring cells using scVI denoised expression...")
    adata_denoised = score_isg_signature(adata_denoised, logger)
    
    # Compare raw vs denoised ISG scores
    if "ISG_score" not in adata.obs.columns:
        adata = score_isg_signature(adata, logger)
    
    adata.obs["ISG_score_raw"] = adata.obs["ISG_score"]
    adata.obs["ISG_score_denoised"] = adata_denoised.obs["ISG_score"]
    
    # Plot comparison
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    sns.scatterplot(
        data=adata.obs, x="ISG_score_raw", y="ISG_score_denoised",
        hue=config.DISEASE_KEY, alpha=0.5, ax=axes[0]
    )
    axes[0].set_title("ISG Score: Raw vs scVI Denoised")
    axes[0].set_xlabel("Raw ISG Score")
    axes[0].set_ylabel("scVI Denoised ISG Score")
    
    # Boxplot of denoised ISG score by disease status in monocytes
    mono_mask = adata.obs[config.CELLTYPE_COARSE_KEY].isin(config.MONOCYTE_TYPES)
    sns.boxplot(
        data=adata.obs[mono_mask], x=config.DISEASE_KEY, y="ISG_score_denoised", ax=axes[1]
    )
    axes[1].set_title("Denoised ISG Score by Disease (Monocytes)")
    
    fig.tight_layout()
    save_figure(fig, "06_isg_score_raw_vs_denoised", logger)
    
    # Correlation between raw and denoised score
    corr = adata.obs["ISG_score_raw"].corr(adata.obs["ISG_score_denoised"])
    logger.info("Pearson correlation between raw and denoised ISG scores: %.3f", corr)
    
    # Save a comparison table
    comparison_df = adata.obs[[config.INDIVIDUAL_KEY, config.DISEASE_KEY, "ISG_score_raw", "ISG_score_denoised"]]
    save_table(comparison_df, "06_isg_score_raw_vs_denoised_comparison", logger)


if __name__ == "__main__":
    main()
