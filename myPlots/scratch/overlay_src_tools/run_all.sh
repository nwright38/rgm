#!/bin/bash
#
# run_all.sh
#
# Run all four SRC analysis macros with hardcoded file/label pairs.
# Edit the FILE/LABEL pairs below, then run: ./run_all.sh
#
# Common cuts and parameters apply to all macros.

set -e

################################################################################
# COMMON CONFIGURATION (applies to all macros)
################################################################################

OUTPUT_DIR="pdf"
TREE_NAME="srcTree"
EPP_CUT="pCM > 0"
BASE_CUT="pCM > 0 && pMiss < 1. && recP < 1."
APP_CUT=""
NORMALIZE_TO_UNITY="true"
INCLUDE_FD_FD="false"
Q2_EDGES="1.5,1.80,2.10,2.40,3.00,5.0"
AUTO_REBIN="true"
MAX_REBIN_FACTOR="2"

################################################################################
# HARDCODED FILE/LABEL PAIRS — Edit these, then run ./run_all.sh
################################################################################

# DATA_FILE and SIM_FILE pairs (comment/uncomment as needed)

# Carbon-12
DATA_FILE="~/data/RGM_DATA/c12_src_skim.root"
SIM_FILE="~/data/RGM_DATA/c12_sim_skim.root"
DATA_LABEL="C12 Data"
SIM_LABEL="C12 Sim"
OUTPUT_PREFIX="c12"  # prefix for output filenames (e.g., c12, he4, etc)

# Helium-4
# DATA_FILE="~/data/RGM_DATA/he4_src_skim.root"
# SIM_FILE="~/data/RGM_DATA/he4_sim_skim_100MeV_allD.root"
# DATA_LABEL="He4 Data"
# SIM_LABEL="He4 Sim"
# OUTPUT_PREFIX="he4"  # prefix for output filenames (e.g., c12, he4, etc)

################################################################################
# END CONFIGURATION
################################################################################

# expand tilde
DATA_FILE="${DATA_FILE/#\~/$HOME}"
SIM_FILE="${SIM_FILE/#\~/$HOME}"

# validate
if [[ ! -f "$DATA_FILE" ]]; then
  echo "Error: Data file not found: $DATA_FILE" >&2
  exit 1
fi
if [[ ! -f "$SIM_FILE" ]]; then
  echo "Error: Sim file not found: $SIM_FILE" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "=================================================="
echo "SRC Analysis Running on:"
echo "  Data: $DATA_LABEL ($DATA_FILE)"
echo "  Sim:  $SIM_LABEL ($SIM_FILE)"
echo "  EPP cut: $EPP_CUT"
echo "  Base cut: $BASE_CUT"
echo "=================================================="
echo ""

# Uncomment individual lines to run specific macros, or remove comments to run all.

# Data single-file analysis
echo "[1/7] overlay_data_by_detector (data by detector)..."
root -l -b -q "overlay_data_by_detector.C(\"$DATA_FILE\",\"$TREE_NAME\",\"${OUTPUT_DIR}/${OUTPUT_PREFIX}_data_overlay_data_by_detector.pdf\",$NORMALIZE_TO_UNITY,\"$EPP_CUT\",\"$BASE_CUT\",\"\",\"$APP_CUT\",$INCLUDE_FD_FD)"
echo ""

# Sim single-file analysis
echo "[2/7] overlay_data_by_detector (sim by detector)..."
root -l -b -q "overlay_data_by_detector.C(\"$SIM_FILE\",\"$TREE_NAME\",\"${OUTPUT_DIR}/${OUTPUT_PREFIX}_sim_overlay_data_by_detector.pdf\",$NORMALIZE_TO_UNITY,\"$EPP_CUT\",\"$BASE_CUT\",\"(weight_epp)*(weight_epp<200)\",\"$APP_CUT\",$INCLUDE_FD_FD)"
echo ""

# Data + Sim overlay
echo "[3/7] overlay_default_multi (data vs sim overlay)..."
FILES_CSV="$DATA_FILE,$SIM_FILE"
WEIGHTS_CSV="1.0,(weight_epp)*(weight_epp<200)"
LABELS_CSV="$DATA_LABEL,$SIM_LABEL"
root -l -b -q "overlay_default_multi.C(\"$FILES_CSV\",\"$TREE_NAME\",\"${OUTPUT_DIR}/${OUTPUT_PREFIX}_data_sim_overlay_default_multi.pdf\",$NORMALIZE_TO_UNITY,\"$EPP_CUT\",\"$BASE_CUT\",\"$WEIGHTS_CSV\",\"$LABELS_CSV\",\"$APP_CUT\",$INCLUDE_FD_FD)"
echo ""

# Data 2D
echo "[4/7] h2d_src_fdcd (data 2D distributions)..."
root -l -b -q "h2d_src_fdcd.C(\"$DATA_FILE\",\"$TREE_NAME\",\"${OUTPUT_DIR}/${OUTPUT_PREFIX}_data_h2d_src_fdcd.pdf\",\"$EPP_CUT\",\"$BASE_CUT\",\"\",\"$APP_CUT\",$INCLUDE_FD_FD)"
echo ""

# Sim 2D
echo "[5/7] h2d_src_fdcd (sim 2D distributions)..."
root -l -b -q "h2d_src_fdcd.C(\"$SIM_FILE\",\"$TREE_NAME\",\"${OUTPUT_DIR}/${OUTPUT_PREFIX}_sim_h2d_src_fdcd.pdf\",\"$EPP_CUT\",\"$BASE_CUT\",\"(weight_epp)*(weight_epp<200)\",\"$APP_CUT\",$INCLUDE_FD_FD)"
echo ""

# Data Q2 bins
echo "[6/7] overlay_q2_by_detector (data Q^2 bins)..."
root -l -b -q "overlay_q2_by_detector.C(\"$DATA_FILE\",\"$TREE_NAME\",\"${OUTPUT_DIR}/${OUTPUT_PREFIX}_data_overlay_q2_by_detector.pdf\",$NORMALIZE_TO_UNITY,\"$EPP_CUT\",\"$BASE_CUT\",\"\",\"$APP_CUT\",\"$Q2_EDGES\",$AUTO_REBIN,$MAX_REBIN_FACTOR,$INCLUDE_FD_FD)"
echo ""

# Sim Q2 bins
echo "[7/7] overlay_q2_by_detector (sim Q^2 bins)..."
root -l -b -q "overlay_q2_by_detector.C(\"$SIM_FILE\",\"$TREE_NAME\",\"${OUTPUT_DIR}/${OUTPUT_PREFIX}_sim_overlay_q2_by_detector.pdf\",$NORMALIZE_TO_UNITY,\"$EPP_CUT\",\"$BASE_CUT\",\"(weight_epp)*(weight_epp<200)\",\"$APP_CUT\",\"$Q2_EDGES\",$AUTO_REBIN,$MAX_REBIN_FACTOR,$INCLUDE_FD_FD)"
echo ""

echo "=================================================="
echo "Done! Check ${OUTPUT_DIR}/ for output PDFs."
echo "=================================================="
