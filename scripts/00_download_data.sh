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
FILE_NAME="GSE174188_CLUES1_adjusted.h5ad.gz"
# NCBI GEO HTTP download URL (used as primary source because the FTP suppl directory is omitted for this accession)
DOWNLOAD_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=${GEO_ACCESSION}&format=file&file=${FILE_NAME}"

DEST_GZ="$RAW_DIR/$FILE_NAME"
DEST_H5AD="$RAW_DIR/${FILE_NAME%.gz}"

if [[ -f "$DEST_H5AD" ]]; then
  echo "Already present: $DEST_H5AD (delete it if you want to re-download)"
  exit 0
fi

echo "Downloading $GEO_ACCESSION supplementary file from GEO..."
echo "  $DOWNLOAD_URL"
curl -L --fail --retry 3 --retry-delay 5 -o "$DEST_GZ" "$DOWNLOAD_URL"

echo "Decompressing..."
gunzip -k "$DEST_GZ"   # -k keeps the .gz in case you want to re-verify later
rm -f "$DEST_GZ"

echo "Done: $DEST_H5AD"
echo ""
echo "If this curl command fails (GEO occasionally reorganizes FTP paths), grab"
echo "the file manually from the GEO Series page instead:"
echo "  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188"
echo "-> Supplementary file: GSE174188_CLUES1_adjusted.h5ad.gz"
echo "and place the decompressed .h5ad at: $DEST_H5AD"
