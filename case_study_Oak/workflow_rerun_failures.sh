#!/bin/bash
# =====================================================================
# Complete workflow to identify and rerun failed traits
# =====================================================================

set -euo pipefail

CASE_DIR="/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/case_study"

cd "$CASE_DIR"


echo "MONET Convergence Failure Workflow"
echo ""

# Step 1: Identify failed traits
echo "Step 1: Identifying traits with convergence issues..."
echo "        (Rhat > 1.05 or ESS < 300)"
echo ""

module load r-light
Rscript identify_convergence_failures.r

echo ""

# Check if any failures were found
if [ ! -f "results_logAV_env/traits_needing_rerun.csv" ]; then
    echo "ERROR: traits_needing_rerun.csv not created."
    echo "Check identify_convergence_failures.r output above."
    exit 1
fi

n_failures=$(tail -n +2 "results_logAV_env/traits_needing_rerun.csv" | wc -l)

if [ "$n_failures" -eq 0 ]; then
    echo "✓ No convergence failures found! All traits look good."
    exit 0
fi

echo ""
echo "========================================="
echo "Step 2: Submitting rerun jobs"
echo "========================================="
echo ""
echo "Found $n_failures trait(s) needing rerun."
echo ""
echo "Submitting SLURM array job with STRICT settings:"
echo "  - warmup=4000, iter=12000, thin=4"
echo "  - adapt_delta=0.99, max_treedepth=15"
echo "  - 72 hour time limit"
echo "  - 32GB memory"
echo ""

sbatch job_rerun_failed_traits.sh

echo ""
echo "Job submitted. Monitor progress with:"
echo "    squeue -u $USER"
echo "    tail -f environment_logs/rerun_*.stdout"
