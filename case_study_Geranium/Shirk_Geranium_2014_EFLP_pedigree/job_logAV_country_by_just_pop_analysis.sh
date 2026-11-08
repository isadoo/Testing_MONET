#!/bin/bash
#SBATCH --mail-user ijeronim@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --job-name=monet_country_by_just_pop
#SBATCH --output=country_by_just_pop_logs/country_by_just_pop_%A_%a.stdout
#SBATCH --error=country_by_just_pop_logs/country_by_just_pop_%A_%a.stderr
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=30:00:00
#SBATCH --account jgoudet_pop_fst
#SBATCH --array=0-4


set -euo pipefail

mkdir -p country_by_just_pop_logs

export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

EXP1_TRAITS=(boltA boltB boltC)
EXP2_TRAITS=(totC totN)

PAIRS=()
for t in "${EXP1_TRAITS[@]}"; do PAIRS+=("exp1:$t"); done
for t in "${EXP2_TRAITS[@]}"; do PAIRS+=("exp2:$t"); done

idx=${SLURM_ARRAY_TASK_ID:-0}
if [ "$idx" -ge "${#PAIRS[@]}" ]; then
  echo "SLURM_ARRAY_TASK_ID=$idx has no matching entry in PAIRS (size ${#PAIRS[@]})." >&2
  exit 1
fi
pair="${PAIRS[$idx]}"
export MONET_EXPERIMENT="${pair%%:*}"
export MONET_TRAIT="${pair##*:}"

CASE_DIR="/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/case_study_Geranium/Shirk_Geranium_2014_EFLP_pedigree"

export MONET_BASEDIR="$CASE_DIR"
export MONET_OUTDIR="$CASE_DIR/results_logAV_country_by_just_pop"

echo "Task $idx: experiment=$MONET_EXPERIMENT trait=$MONET_TRAIT"

module load r-light

cd "$CASE_DIR"
Rscript logAV_country_by_just_pop_analysis.r
