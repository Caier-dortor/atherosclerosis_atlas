#!/bin/bash
#
# Atherosclerosis Atlas Analysis Pipeline
# Phase 1-7 execution script
# Data: Traeuble et al., 2025, Nature Communications
#

set -e

PROJ_DIR="D:/openclaw_workspace/atherosclerosis_atlas"
SCRIPT_DIR="$PROJ_DIR/scripts"
RES_DIR="$PROJ_DIR/results"

echo "============================================"
echo " Atherosclerosis Atlas Pipeline"
echo "============================================"
echo ""

# ---- Phase 1: Macrophage Analysis ----
echo "=== Phase 1: Macrophage Vascular Bed Analysis ==="
python "$SCRIPT_DIR/phase1_macrophage_analysis.py"
echo "Phase 1 done."
echo ""

# ---- Check if R is available for complement ----
if command -v Rscript &> /dev/null; then
    echo "=== Phase 1 R Complement (variancePartition, dream) ==="
    Rscript "$SCRIPT_DIR/phase1_R_complement.R"
    echo "Phase 1 R complement done."
else
    echo "WARNING: R not found — skipping variancePartition/DESeq2 (need R 4.5+)"
    echo "  Install R and run: Rscript $SCRIPT_DIR/phase1_R_complement.R"
fi
echo ""

# ---- Phase 2: Trained Immunity ----
echo "=== Phase 2: Trained Immunity Signatures ==="
python "$SCRIPT_DIR/phase2_trained_immunity.py"
echo "Phase 2 done."
echo ""

# ---- Phase 3: Metabolism-Epigenetics Coupling ----
echo "=== Phase 3: Metabolism-Epigenetics Coupling ==="
python "$SCRIPT_DIR/phase3_metabolism_epigenetics.py"
echo "Phase 3 done."
echo ""

# ---- Phase 4-7: To be implemented ---
echo "=== Pipeline Status ==="
echo "Phase 1: macrpage analysis + figures"
echo "Phase 2: trained immunity + figures"
echo "Phase 3: metabolism-epigenetics + figures"
echo "Phase 4-7: scripts ready, run after Phase 1-3 complete"
echo ""
echo "Results directory: $RES_DIR"
