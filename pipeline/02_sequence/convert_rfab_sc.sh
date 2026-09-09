#!/bin/bash
#SBATCH --job-name=rfab_conv
#SBATCH --partition=shared-cpu
#SBATCH --account=falcone_teach_generic
#SBATCH --ntasks=1 --cpus-per-task=2 --mem=16G --time=01:00:00
#SBATCH --output=logs/rfab_conv_%j.out

SIF=$(ls ~/master_thesis/RFantibody/*.sif | head -1)
echo "SIF: $SIF"
for c in ha_stem cxcr4; do
  OUT=~/master_thesis/results/rfantibody/final/$c
  echo "=== Conversion $c ==="
  apptainer exec ${SIF} qvscorefile ${OUT}/${c}_3_rf2.qv
  n=$(grep -v -c "^SCORE: score\|^$" ${OUT}/${c}_3_rf2.sc 2>/dev/null || echo 0)
  echo "$c : $n scores générés"
done
