#!/bin/bash
set -euo pipefail

# Draft Reference Pipeline — reproducible subset (steps 1-4 and 8)
# Assumes this script is run from the repository root.

cd "$(dirname "$0")"

SCRIPTS=scripts

echo "=========================================="
echo "Step 1: MS harmonization"
echo "=========================================="
python "$SCRIPTS/harmonize_massspec.py"

echo ""
echo "=========================================="
echo "Step 2: Generate parameter ranges"
echo "=========================================="
python "$SCRIPTS/generate_parameter_ranges.py"
python "$SCRIPTS/generate_tRNA_parameter_ranges.py"

echo ""
echo "=========================================="
echo "Step 3: Grid search"
echo "=========================================="
python "$SCRIPTS/run_rrna_grid_from_ranges.py"
python "$SCRIPTS/run_3way_nm_from_ranges.py"
python "$SCRIPTS/run_3way_m5C_from_ranges.py"
python "$SCRIPTS/run_3way_m6A_from_ranges.py"
python "$SCRIPTS/run_polya_grid_fast.py"
python "$SCRIPTS/run_trna_grid_from_ranges.py"

echo ""
echo "=========================================="
echo "Step 4: Filtering"
echo "=========================================="
python "$SCRIPTS/run_filtering.py"
python "$SCRIPTS/run_tRNA_filtering.py"

echo ""
echo "=========================================="
echo "Step 8: Tiered lists"
echo "=========================================="
python "$SCRIPTS/create_tiered_mod_lists.py"
python "$SCRIPTS/create_tRNA_tiered_table.py"

echo ""
echo "=========================================="
echo "Pipeline complete (steps 1-4 and 8)."
echo "=========================================="
