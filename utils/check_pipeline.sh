#!/bin/bash
T="mdm2 pdl1 spike tnfa egfr insulinr vegf stat3 kit glp1r pdgfrb fgfr2 trka ha_head ha_stem il7ra clathrin cd28 ccr5 cxcr4"

echo "############ ETAT PIPELINE $(date +%H:%M) ############"
echo ""
echo "=== JOBS QUI TOURNENT ==="
squeue -u syrvet4 -h -o "%j" | sed 's/_[0-9]*$//' | sort | uniq -c
echo ""

printf "%-9s | %-11s | %-9s | %-9s | %-9s | %-9s\n" "CIBLE" "BENNETT" "PPDIFF-of3" "BC-of3" "PPD-pyro" "BC-pyro"
echo "-----------------------------------------------------------------------------"
for c in $T; do
  # Bennett : mpnn / af2
  bm=$(ls ~/master_thesis/results/bennett/official/$c/mpnn_out/*.pdb 2>/dev/null | wc -l)
  ba=$(ls ~/master_thesis/results/bennett/official/$c/af2_out/*.pdb 2>/dev/null | wc -l)
  # OF3
  po=$(ls -d ~/master_thesis/shared/OpenFold3/outputs/ppdiff/$c/*/ 2>/dev/null | wc -l)
  bo=$(ls -d ~/master_thesis/shared/OpenFold3/outputs/bindcraft/$c/*/ 2>/dev/null | wc -l)
  # PyRosetta (lignes du CSV - 1)
  pp=$([ -f ~/master_thesis/results/ranking/ppdiff_pyrosetta/${c}_pyrosetta.csv ] && echo $(($(wc -l < ~/master_thesis/results/ranking/ppdiff_pyrosetta/${c}_pyrosetta.csv)-1)) || echo "-")
  bp=$([ -f ~/master_thesis/results/ranking/bindcraft_pyrosetta/${c}_pyrosetta.csv ] && echo $(($(wc -l < ~/master_thesis/results/ranking/bindcraft_pyrosetta/${c}_pyrosetta.csv)-1)) || echo "-")
  printf "%-9s | %4s/%-4s | %4s     | %4s    | %4s     | %4s\n" "$c" "$bm" "$ba" "$po" "$bo" "$pp" "$bp"
done

echo ""
echo "=== BOLTZGEN (dfold/CSV) ==="
for c in $T; do
  df=$(ls ~/master_thesis/results/boltzgen/final/$c/intermediate_designs_inverse_folded/fold_out_design_npz/*.npz 2>/dev/null | wc -l)
  csv=$([ -f ~/master_thesis/results/boltzgen/final/$c/final_ranked_designs/all_designs_metrics.csv ] && echo "CSV" || echo "-")
  printf "%-9s dfold=%5s %s\n" "$c" "$df" "$csv"
done

echo ""
echo "=== TOTAUX ==="
tm=0; ta=0
for c in $T; do
  tm=$((tm+$(ls ~/master_thesis/results/bennett/official/$c/mpnn_out/*.pdb 2>/dev/null | wc -l)))
  ta=$((ta+$(ls ~/master_thesis/results/bennett/official/$c/af2_out/*.pdb 2>/dev/null | wc -l)))
done
echo "Bennett mpnn: $tm/20000 | af2: $ta/20000"
echo "BoltzGen CSV: $(ls ~/master_thesis/results/boltzgen/final/*/final_ranked_designs/all_designs_metrics.csv 2>/dev/null | wc -l)/20"
echo "PyRosetta PPDiff CSV: $(ls ~/master_thesis/results/ranking/ppdiff_pyrosetta/*_pyrosetta.csv 2>/dev/null | wc -l)/20"
echo "PyRosetta BindCraft CSV: $(ls ~/master_thesis/results/ranking/bindcraft_pyrosetta/*_pyrosetta.csv 2>/dev/null | wc -l)/15"
