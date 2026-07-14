"""
Shared constants for the lupus batch-integration pipeline.

Nothing in here should require scanpy/anndata to import — keep this module
lightweight so tests can import it without the full environment installed.
"""
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = REPO_ROOT / "data" / "raw"
DATA_PROCESSED = REPO_ROOT / "data" / "processed"
RESULTS_FIGURES = REPO_ROOT / "results" / "figures"
RESULTS_TABLES = REPO_ROOT / "results" / "tables"

RAW_H5AD = DATA_RAW / "GSE174188_CLUES1_adjusted.h5ad"

CKPT_DEV_SUBSET = DATA_PROCESSED / "lupus_dev_subset.h5ad"
CKPT_QC = DATA_PROCESSED / "lupus_qc.h5ad"
CKPT_BASELINE = DATA_PROCESSED / "lupus_baseline_embedding.h5ad"
CKPT_HARMONY = DATA_PROCESSED / "lupus_harmony.h5ad"
CKPT_BBKNN = DATA_PROCESSED / "lupus_bbknn.h5ad"
CKPT_SCVI = DATA_PROCESSED / "lupus_scvi.h5ad"

# ---------------------------------------------------------------------------
# Metadata field names (see README "Key metadata fields" table)
# ---------------------------------------------------------------------------
BATCH_KEY = "batch_cov"                # primary technical batch: mux-seq pool
PROCESSING_COHORT_KEY = "Processing_Cohort"  # coarser processing-time batch
INDIVIDUAL_KEY = "ind_cov"
SAMPLE_KEY = "ind_cov_batch_cov"       # individual x batch, unique per sample
DISEASE_KEY = "SLE_status"             # "SLE" / "Healthy"
STATUS_KEY = "Status"
ANCESTRY_KEY = "pop_cov"
CELLTYPE_COARSE_KEY = "cg_cov"         # published ground-truth label, 11 types
CELLTYPE_FINE_KEY = "ct_cov"

# The 11 published coarse cell types (Perez et al. 2022, Fig. 1A)
CELLTYPES_COARSE = [
    "B", "T4", "T8", "NK", "cM", "ncM", "cDC", "pDC", "PB", "Prolif", "Progen",
]
MONOCYTE_TYPES = ["cM", "ncM"]

# ---------------------------------------------------------------------------
# QC thresholds
# ---------------------------------------------------------------------------
# The GEO object is already lightly QC'd by the original authors (they only
# deposited the ~2,000-gene HVG subset with cells that passed their filters).
# These thresholds are a second, independent pass — deliberately conservative
# so we don't discard cells the original pipeline already vetted. Documented
# per the workflow guide's convention of stating + justifying thresholds
# rather than treating them as a universal default.
MIN_GENES_PER_CELL = 200
MIN_CELLS_PER_GENE = 3
MAX_PCT_MT = 20.0  # adult PBMC convention; this panel has few/no MT genes but
                    # we compute it if present for completeness

# ---------------------------------------------------------------------------
# Embedding / clustering parameters
# ---------------------------------------------------------------------------
N_PCS = 30
N_NEIGHBORS = 15
LEIDEN_RESOLUTION = 1.0
RANDOM_STATE = 0

# ---------------------------------------------------------------------------
# Interferon-stimulated gene (ISG) panel
# ---------------------------------------------------------------------------
# A standard, widely used core type-I IFN signature panel (genes recur across
# many independent IFN-signature studies in lupus and other autoimmune
# disease literature, e.g. Baechler et al. 2003 PNAS; Feng et al. 2006
# Arthritis Rheum; Kirou et al. 2004 Arthritis Rheum). This is a generic,
# literature-standard panel, not a reproduction of any single paper's exact
# gene list, and functions as a starting point — swap in
# MSigDB HALLMARK_INTERFERON_ALPHA_RESPONSE for a fully unbiased alternative
# (that comparison is exactly what the GSEA step in script 06 does).
CORE_ISG_PANEL = [
    "ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT3",
    "MX1", "MX2", "OAS1", "OAS2", "OAS3", "RSAD2", "STAT1", "IRF7",
    "USP18", "LY6E", "XAF1", "HERC5", "EIF2AK2",
]

# Downsampling defaults for the local development subset
DEV_N_CELLS_DEFAULT = 120_000
DEV_SEED = 0
