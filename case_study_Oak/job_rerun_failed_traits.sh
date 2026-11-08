#!/bin/bash
#SBATCH --mail-user ijeronim@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --job-name=monet_rerun
#SBATCH --output=environment_logs/rerun_%A_%a.stdout
#SBATCH --error=environment_logs/rerun_%A_%a.stderr
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --account jgoudet_pop_fst
#SBATCH --array=0-2

# ---------------------------------------------------------------
# Rerun problematic traits with STRICT MCMC settings
# 
# First run identify_convergence_failures.r to generate
# the traits_needing_rerun.csv file, then submit this job.
# ---------------------------------------------------------------

set -euo pipefail

mkdir -p environment_logs

SPECIES_ARR=(faginea faginea lusitanica)
GARDEN_ARR=(GH      OC      GH)

idx=${SLURM_ARRAY_TASK_ID:-0}
export MONET_SPECIES="${SPECIES_ARR[$idx]}"
export MONET_GARDEN="${GARDEN_ARR[$idx]}"

CASE_DIR="/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/case_study"
TOOLS_DIR="/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/tools"

export MONET_BASEDIR="$CASE_DIR"
export MONET_TOOLSDIR="$TOOLS_DIR"
export MONET_OUTDIR="$CASE_DIR/results_logAV_env"

echo "Rerun task $idx: species=$MONET_SPECIES garden=$MONET_GARDEN"
echo "STRICT MCMC: warmup=4000, iter=12000, thin=4, adapt_delta=0.99"

module load r-light

cd "$CASE_DIR"
Rscript logAV_rerun_failed_traits.r
