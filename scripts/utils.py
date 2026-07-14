"""
Shared helpers used across pipeline scripts: logging, checkpointing, QC,
and small plotting utilities. Kept separate from config.py so that config
stays importable without scanpy installed.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Iterable

import config


def get_logger(name: str) -> logging.Logger:
    """Consistent logger so every script's stdout is uniformly formatted and
    timestamped — makes it easy to tell how long each pipeline stage took
    when rerunning on the full 1.26M-cell object."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(
            logging.Formatter("[%(asctime)s] %(name)s: %(message)s", "%H:%M:%S")
        )
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger


def require_file(path: Path, hint: str = "") -> None:
    """Fail fast with a useful message instead of a bare FileNotFoundError
    several stack frames deep inside scanpy's h5py reader."""
    if not path.exists():
        msg = f"Expected input file not found: {path}"
        if hint:
            msg += f"\n  -> {hint}"
        raise FileNotFoundError(msg)


def resolve_gene_panel(adata, panel: Iterable[str], logger: logging.Logger | None = None):
    """This dataset's expression matrix only contains the ~1,999-gene HVG
    subset the original authors used (see README "Data note"), so any curated
    gene panel needs to be intersected with what's actually present before
    scoring. Returns the subset of `panel` found in `adata.var_names` and logs
    which genes were dropped, rather than silently scoring on fewer genes than
    the caller thinks they asked for.
    """
    panel = list(panel)
    present = [g for g in panel if g in adata.var_names]
    missing = sorted(set(panel) - set(present))
    if logger and missing:
        logger.info(
            "Gene panel: %d/%d genes present in this object; missing: %s",
            len(present), len(panel), ", ".join(missing),
        )
    if len(present) < 3:
        raise ValueError(
            f"Only {len(present)} of the requested panel genes are present in "
            "this object — too few to compute a meaningful signature score."
        )
    return present


def qc_metrics(adata, logger: logging.Logger | None = None):
    """Compute standard QC metrics. Mitochondrial genes are flagged if present,
    but note this object is pre-restricted to ~1,999 HVGs, so %MT here is not
    directly comparable to a %MT computed on the full transcriptome — it's
    included for completeness/consistency with standard practice, not as a
    primary filtering criterion for this particular object.
    """
    import scanpy as sc

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    n_mt = int(adata.var["mt"].sum())
    if logger:
        logger.info("Mitochondrial genes present in panel: %d", n_mt)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"] if n_mt > 0 else [], percent_top=None,
        log1p=False, inplace=True,
    )
    return adata


def save_checkpoint(adata, path: Path, logger: logging.Logger | None = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(path)
    if logger:
        logger.info("Wrote checkpoint: %s (%d cells x %d genes)",
                     path, adata.n_obs, adata.n_vars)


def setup_plot_style():
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="notebook")
    plt.rcParams["figure.dpi"] = 120
    plt.rcParams["savefig.bbox"] = "tight"


def save_figure(fig, name: str, logger: logging.Logger | None = None):
    config.RESULTS_FIGURES.mkdir(parents=True, exist_ok=True)
    out_path = config.RESULTS_FIGURES / f"{name}.png"
    fig.savefig(out_path)
    if logger:
        logger.info("Saved figure: %s", out_path)
    return out_path


def save_table(df, name: str, logger: logging.Logger | None = None):
    config.RESULTS_TABLES.mkdir(parents=True, exist_ok=True)
    out_path = config.RESULTS_TABLES / f"{name}.csv"
    df.to_csv(out_path)
    if logger:
        logger.info("Saved table: %s", out_path)
    return out_path
