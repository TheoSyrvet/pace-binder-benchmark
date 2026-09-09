#!/bin/bash
#SBATCH --job-name=pyro_bc
#SBATCH --partition=shared-cpu
#SBATCH --account=falcone_teach_generic
#SBATCH --ntasks=1 --cpus-per-task=4 --mem=16G --time=12:00:00
#SBATCH --array=1-20
#SBATCH --output=logs/pyro_bc_%a_%j.out

source $(conda info --base)/etc/profile.d/conda.sh
conda activate af2_binder_design
export PYTHONNOUSERSITE=1

TARGETS=("" "mdm2" "pdl1" "spike" "tnfa" "egfr" "insulinr" "vegf" "stat3" "kit" "glp1r" "pdgfrb" "fgfr2" "trka" "ha_head" "ha_stem" "il7ra" "clathrin" "cd28" "ccr5" "cxcr4")
TARGET=${TARGETS[$SLURM_ARRAY_TASK_ID]}

ACC=~/master_thesis/results/bindcraft/final/${TARGET}/Accepted
AGG=~/master_thesis/results/ranking/bindcraft_pyrosetta/agg/${TARGET}_agg.csv
OUT=~/master_thesis/results/ranking/bindcraft_pyrosetta/${TARGET}_pyrosetta.csv

[ ! -f "$AGG" ] && { echo "pas d'agg pour $TARGET, skip"; exit 0; }
mkdir -p ~/master_thesis/shared/scripts/logs

echo "=== PyRosetta BindCraft ${TARGET} (binder=B, target=A) ==="
python3 ~/master_thesis/shared/scripts/pyrosetta_bindcraft.py \
    --colabfold_dir "$ACC" \
    --agg_csv "$AGG" \
    --target "$TARGET" \
    --binder_chain B \
    --target_chain A \
    --tool_filter all \
    --output "$OUT"
echo "=== DONE ${TARGET} ==="
