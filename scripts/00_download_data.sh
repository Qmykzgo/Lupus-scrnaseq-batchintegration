#!/usr/bin/env bash
# Downloads the processed AnnData object for GSE174188 (Perez et al. 2022,
# Science) directly from the NCBI GEO FTP mirror. No login required — this is
# open-access data; only the linked genotype data (dbGaP phs002812) is
# access-controlled, and this project never touches that.
#
# Output: data/raw/GSE174188_CLUES1_adjusted.h5ad  (~2-4 GB uncompressed)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
RAW_DIR="$REPO_ROOT/data/raw"
mkdir -p "$RAW_DIR"

GEO_ACCESSION="GSE174188"
FILE_NAME="GSE174188_CLUES1_adjusted.h5ad"
# CELLxGENE direct public download URL (since original GEO FTP/HTTP supplementary downloads are restricted/deleted)
DOWNLOAD_URL="https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad"

DEST_H5AD="$RAW_DIR/$FILE_NAME"

if [[ -f "$DEST_H5AD" ]]; then
  echo "Already present: $DEST_H5AD (delete it if you want to re-download)"
  exit 0
fi

echo "Downloading $GEO_ACCESSION processed dataset from CELLxGENE..."
echo "  $DOWNLOAD_URL"
curl -L --fail --retry 3 --retry-delay 5 -o "$DEST_H5AD" "$DOWNLOAD_URL"

echo "Done: $DEST_H5AD"
echo ""
echo "If this curl command fails (GEO occasionally reorganizes FTP paths), grab"
echo "the file manually from the GEO Series page instead:"
echo "  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188"
echo "-> Supplementary file: GSE174188_CLUES1_adjusted.h5ad.gz"
echo "and place the decompressed .h5ad at: $DEST_H5AD"
