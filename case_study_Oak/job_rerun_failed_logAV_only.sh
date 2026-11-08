#!/bin/bash
#SBATCH --mail-user ijeronim@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --job-name=monet_rerun_logav
#SBATCH --output=rerun_logAV/rerun_%A_%a.stdout
#SBATCH --error=rerun_logAV/rerun_%A_%a.stderr
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --account jgoudet_pop_fst
#SBATCH --array=0-0

mkdir -p rerun_logAV

SPECIES_ARR=(lusitanica)
#GARDEN_ARR=(GH      OC      GH)

idx=${SLURM_ARRAY_TASK_ID:-0}
export MONET_SPECIES="${SPECIES_ARR[$idx]}"
#export MONET_GARDEN="${GARDEN_ARR[$idx]}"
MONET_GARDEN="GH"

CASE_DIR="/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/case_study"
TOOLS_DIR="/work/FAC/FBM/DEE/jgoudet/default/isaChapter2/isaChapter2/Testing_simulated_data/Testing_MONET/tools"

export MONET_BASEDIR="$CASE_DIR"
export MONET_TOOLSDIR="$TOOLS_DIR"
export MONET_OUTDIR="$CASE_DIR/results_logAV_only"

echo "Task $idx: species=$MONET_SPECIES garden=$MONET_GARDEN"

module load r-light

cd "$CASE_DIR"
Rscript logAV_only_rerun_failed_traits.r