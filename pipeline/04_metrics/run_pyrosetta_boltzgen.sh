#!/bin/bash
#SBATCH --job-name=pyr_bg
#SBATCH --partition=shared-cpu
#SBATCH --account=falcone_teach_generic
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=1-20
#SBATCH --output=logs/pyr_bg_%a_%j.out

TARGETS=("" mdm2 pdl1 spike tnfa egfr insulinr vegf stat3 kit glp1r pdgfrb fgfr2 trka ha_head ha_stem il7ra clathrin cd28 ccr5 cxcr4)
TARGET=${TARGETS[$SLURM_ARRAY_TASK_ID]}

# kit abandonne
if [ "$TARGET" = "kit" ]; then echo "kit abandonne, skip"; exit 0; fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate BindCraft
echo ">>> PyRosetta BoltzGen : $TARGET"
python3 -u ~/master_thesis/shared/scripts/pyrosetta_boltzgen.py --target $TARGET
echo ">>> $TARGET termine"
