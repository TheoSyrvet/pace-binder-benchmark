#!/bin/bash
#SBATCH --job-name=pyro_ppd
#SBATCH --partition=shared-cpu
#SBATCH --account=falcone_teach_generic
#SBATCH --ntasks=1 --cpus-per-task=4 --mem=16G --time=12:00:00
#SBATCH --array=1-20
#SBATCH --output=logs/pyro_ppd_%a_%j.out

source $(conda info --base)/etc/profile.d/conda.sh
conda activate af2_binder_design
export PYTHONNOUSERSITE=1

# DB PyRosetta : lecture directe depuis conda (pas de copie, evite saturation I/O en parallele)

TARGETS=("" "mdm2" "pdl1" "spike" "tnfa" "egfr" "insulinr" "vegf" "stat3" "kit" "glp1r" "pdgfrb" "fgfr2" "trka" "ha_head" "ha_stem" "il7ra" "clathrin" "cd28" "ccr5" "cxcr4")
TARGET=${TARGETS[$SLURM_ARRAY_TASK_ID]}

CF_DIR=~/master_thesis/results/colabfold/ppdiff/${TARGET}
AGG=~/master_thesis/results/ranking/ppdiff_pyrosetta/agg/${TARGET}_agg.csv
OUT=~/master_thesis/results/ranking/ppdiff_pyrosetta/${TARGET}_pyrosetta.csv

mkdir -p ~/master_thesis/shared/scripts/logs

echo "=== PyRosetta PPDiff ${TARGET} (binder=A, target=B) ==="
python3 ~/master_thesis/shared/scripts/pyrosetta_post_colabfold_v2.py \
    --colabfold_dir "$CF_DIR" \
    --agg_csv "$AGG" \
    --target "$TARGET" \
    --binder_chain A \
    --target_chain B \
    --tool_filter all \
    --output "$OUT"
echo "=== DONE ${TARGET} ==="
