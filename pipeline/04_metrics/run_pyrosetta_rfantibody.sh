#!/bin/bash
#SBATCH --job-name=pyr_rfab
#SBATCH --partition=shared-cpu
#SBATCH --account=falcone_teach_generic
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=1-20
#SBATCH --output=logs/pyr_rfab_%a_%j.out

TARGETS=("" mdm2 pdl1 spike tnfa egfr insulinr vegf stat3 kit glp1r pdgfrb fgfr2 trka ha_head ha_stem il7ra clathrin cd28 ccr5 cxcr4)
TARGET=${TARGETS[$SLURM_ARRAY_TASK_ID]}

source $(conda info --base)/etc/profile.d/conda.sh
conda activate BindCraft

echo ">>> PyRosetta RFantibody : $TARGET"
python3 -u ~/master_thesis/shared/scripts/pyrosetta_rfantibody.py --target $TARGET
echo ">>> $TARGET termine"
