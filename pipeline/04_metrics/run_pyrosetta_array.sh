#!/bin/bash
#SBATCH --job-name=pyros_final
#SBATCH --account=falcone_teach_generic
#SBATCH --partition=shared-cpu
#SBATCH --array=1-20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/pyrosetta_%A_%a.out

TARGETS=("" "mdm2" "pdl1" "spike" "tnfa" "egfr" "insulinr" "vegf" "stat3" "kit"
         "glp1r" "pdgfrb" "fgfr2" "trka" "ha_head" "ha_stem" "il7ra" "clathrin"
         "cd28" "ccr5" "cxcr4")

declare -A BINDER_CHAIN
declare -A TARGET_CHAIN
BINDER_CHAIN["mdm2"]="B";     TARGET_CHAIN["mdm2"]="A"
BINDER_CHAIN["pdl1"]="B";     TARGET_CHAIN["pdl1"]="A"
BINDER_CHAIN["spike"]="A";    TARGET_CHAIN["spike"]="B"
BINDER_CHAIN["tnfa"]="C";     TARGET_CHAIN["tnfa"]="A,B"
BINDER_CHAIN["egfr"]="B";     TARGET_CHAIN["egfr"]="A"
BINDER_CHAIN["insulinr"]="B"; TARGET_CHAIN["insulinr"]="A"
BINDER_CHAIN["vegf"]="C";     TARGET_CHAIN["vegf"]="A,B"
BINDER_CHAIN["stat3"]="B";    TARGET_CHAIN["stat3"]="A"
BINDER_CHAIN["kit"]="B";      TARGET_CHAIN["kit"]="A"
BINDER_CHAIN["glp1r"]="A";    TARGET_CHAIN["glp1r"]="B"
BINDER_CHAIN["pdgfrb"]="B";   TARGET_CHAIN["pdgfrb"]="A"
BINDER_CHAIN["fgfr2"]="B";    TARGET_CHAIN["fgfr2"]="A"
BINDER_CHAIN["trka"]="B";     TARGET_CHAIN["trka"]="A"
BINDER_CHAIN["ha_head"]="B";  TARGET_CHAIN["ha_head"]="A"
BINDER_CHAIN["ha_stem"]="C";  TARGET_CHAIN["ha_stem"]="A,B"
BINDER_CHAIN["il7ra"]="A";    TARGET_CHAIN["il7ra"]="B"
BINDER_CHAIN["clathrin"]="B"; TARGET_CHAIN["clathrin"]="A"
BINDER_CHAIN["cd28"]="A";     TARGET_CHAIN["cd28"]="B"
BINDER_CHAIN["ccr5"]="B";     TARGET_CHAIN["ccr5"]="A"
BINDER_CHAIN["cxcr4"]="B";    TARGET_CHAIN["cxcr4"]="A"

TARGET=${TARGETS[$SLURM_ARRAY_TASK_ID]}
BC=${BINDER_CHAIN[$TARGET]}
TC=${TARGET_CHAIN[$TARGET]}

COLABFOLD_DIR=~/master_thesis/results/colabfold/final/${TARGET}
OUTPUT_DIR=~/master_thesis/results/pyrosetta/final
AGG_CSV=${OUTPUT_DIR}/${TARGET}_agg_colabfold.csv
OUTPUT_CSV=${OUTPUT_DIR}/${TARGET}_pyrosetta.csv

mkdir -p ${OUTPUT_DIR}
mkdir -p ~/master_thesis/shared/scripts/logs

PYTHON="${PACE_PYTHON:?definir PACE_PYTHON, ex: ~/.conda/envs/BindCraft/bin/python}"

# Skip si déjà fait
if [ -f "${OUTPUT_CSV}" ]; then
    echo ">>> ${TARGET} déjà scoré — skip"
    exit 0
fi

echo ">>> Étape 1 — Génération CSV agrégé ColabFold pour ${TARGET}"
${PYTHON} ~/master_thesis/shared/scripts/generate_agg_csv.py \
    --target        ${TARGET} \
    --colabfold_dir ${COLABFOLD_DIR} \
    --output        ${AGG_CSV}

if [ ! -f "${AGG_CSV}" ]; then
    echo "ERREUR: CSV agrégé non généré pour ${TARGET}"
    exit 1
fi

echo ">>> Étape 2 — PyRosetta scoring pour ${TARGET} (binder=${BC}, target=${TC})"
${PYTHON} ~/master_thesis/shared/scripts/pyrosetta_post_colabfold_v2.py \
    --colabfold_dir ${COLABFOLD_DIR} \
    --agg_csv       ${AGG_CSV} \
    --target        ${TARGET} \
    --binder_chain  ${BC} \
    --target_chain  ${TC} \
    --output        ${OUTPUT_CSV}

echo ">>> ${TARGET} terminé — output: ${OUTPUT_CSV}"
