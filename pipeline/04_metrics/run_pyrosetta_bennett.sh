#!/bin/bash
#SBATCH --job-name=pyr_bnt
#SBATCH --partition=shared-cpu
#SBATCH --account=falcone_teach_generic
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=1-20
#SBATCH --output=logs/pyr_bnt_%a_%j.out

TARGETS=("" mdm2 pdl1 spike tnfa egfr insulinr vegf stat3 kit glp1r pdgfrb fgfr2 trka ha_head ha_stem il7ra clathrin cd28 ccr5 cxcr4)
TARGET=${TARGETS[$SLURM_ARRAY_TASK_ID]}

# 3 trimeres exclus de Bennett (pas de top.txt)
case "$TARGET" in
  tnfa|vegf|ha_stem) echo "$TARGET exclu de Bennett (trimere), skip"; exit 0 ;;
esac

# skip si pas de top.txt (cible non generee)
TOPFILE=~/master_thesis/results/ranking/bennett_top100/${TARGET}_top.txt
if [ ! -f "$TOPFILE" ]; then echo "$TARGET : pas de top.txt, skip"; exit 0; fi

source $(conda info --base)/etc/profile.d/conda.sh
conda activate BindCraft
echo ">>> PyRosetta Bennett : $TARGET"
python3 -u ~/master_thesis/shared/scripts/pyrosetta_bennett.py --target $TARGET
