#!/usr/bin/env python3
"""
PyRosetta post-ColabFold scorer
Calcule les métriques Rosetta sur les PDB rank_001 de ColabFold
pour être comparable aux métriques BindCraft.

Usage:
    python3 pyrosetta_post_colabfold.py \
        --colabfold_dir /path/to/colabfold/outputs/complex/mdm2 \
        --agg_csv /path/to/aggregate_mdm2.csv \
        --target mdm2 \
        --binder_chain B \
        --target_chain A \
        --output /path/to/pyrosetta_mdm2_final.csv
"""

import os, sys, csv, re, time, argparse
import pyrosetta
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking

def score_af2_complex(pdb_path, binder_chain, target_chain, sfxn, fr):
    """Score un complexe AF2 avec PyRosetta."""

    pose = pose_from_pdb(pdb_path)

    # Identifier les numéros de chaîne dans la pose
    chain_map = {}
    for i in range(1, pose.num_chains()+1):
        chain_id = pose.pdb_info().chain(pose.chain_begin(i))
        chain_map[chain_id] = i
    
    if binder_chain not in chain_map or not all(c in chain_map for c in target_chain.split(",")):
        raise ValueError(f"Chaînes {binder_chain}/{target_chain} non trouvées. "
                        f"Disponibles: {list(chain_map.keys())}")

    # Pack rotamers pour s'assurer que les chaînes latérales sont optimisées
    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())
    PackRotamersMover(sfxn,
        tf.create_task_and_apply_taskoperations(pose)).apply(pose)

    # FastRelax 1 cycle
    fr.apply(pose)

    # Score binder seul (pour Binder_Energy_Score)
    binder_chain_num = chain_map[binder_chain]
    binder_start = pose.chain_begin(binder_chain_num)
    binder_end   = pose.chain_end(binder_chain_num)
    binder_len   = binder_end - binder_start + 1

    # Extraire le binder seul
    from pyrosetta.rosetta.core.pose import Pose
    from pyrosetta.rosetta.protocols.grafting import delete_region
    binder_pose = Pose(pose, binder_start, binder_end)
    binder_energy = sfxn(binder_pose)
    binder_energy_per_res = binder_energy / max(binder_len, 1)

    # Surface hydrophobicité du binder
    hydrophobic_aa = set("VILMFYW")
    binder_seq = "".join(pose.residue(i).name1()
                         for i in range(binder_start, binder_end+1))
    surface_hydrophobicity = sum(1 for aa in binder_seq
                                  if aa in hydrophobic_aa) / max(len(binder_seq), 1)

    # Interface analysis sur le complexe relaxé
    interface_str = f"{binder_chain}_{''.join(target_chain.split(','))}"
    ia = InterfaceAnalyzerMover()
    ia.set_interface(interface_str)
    ia.set_scorefunction(sfxn)
    ia.set_compute_interface_energy(True)
    ia.set_compute_separated_sasa(True)
    ia.set_compute_interface_sc(True)
    ia.set_compute_interface_delta_hbond_unsat(True)
    ia.set_pack_input(True)
    ia.set_pack_separated(True)
    ia.apply(pose)

    dG       = ia.get_interface_dG()
    dSASA    = ia.get_interface_delta_sasa()
    n_res    = ia.get_num_interface_residues()
    unsat_hb = ia.get_interface_delta_hbond_unsat()
    packstat = ia.get_interface_packstat()
    hbond_e  = ia.get_total_Hbond_E()
    complex_e = ia.get_complex_energy()
    sep_e     = ia.get_separated_interface_energy()

    # Shape Complementarity
    scf = ShapeComplementarityFilter()
    scf.jump_id(1)
    sc_val = scf.score(pose)

    # dG/dSASA ratio (normalisation)
    dg_dsasa = dG / max(dSASA, 1)

    return {
        "dG":                    round(dG, 3),
        "dSASA":                 round(dSASA, 1),
        "SC":                    round(sc_val, 4),
        "n_InterfaceResidues":   n_res,
        "n_InterfaceHbonds":     round(hbond_e, 3),
        "n_InterfaceUnsatHbonds":round(unsat_hb, 1),
        "PackStat":              round(packstat, 4),
        "Binder_Energy_Score":   round(binder_energy, 3),
        "Binder_Energy_per_res": round(binder_energy_per_res, 3),
        "Surface_Hydrophobicity":round(surface_hydrophobicity, 4),
        "dG_dSASA_ratio":        round(dg_dsasa, 6),
        "complex_energy":        round(complex_e, 3),
        "binder_len":            binder_len,
        "binder_seq":            binder_seq,
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--colabfold_dir", required=True)
    parser.add_argument("--agg_csv",       required=True)
    parser.add_argument("--target",        required=True)
    parser.add_argument("--binder_chain",  default="B")
    parser.add_argument("--target_chain",  default="A")
    parser.add_argument("--output",        required=True)
    parser.add_argument("--tool_filter",   default="rfdiff_mpnn",
                        help="Filtrer par outil dans le CSV agrégé (rfdiff_mpnn ou all)")
    args = parser.parse_args()

    pyrosetta.init("-mute all -ex1 -ex2aro")
    sfxn = get_score_function(True)
    fr   = FastRelax(sfxn, 1)

    # Charger le CSV agrégé
    agg_data = {}
    with open(args.agg_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if args.tool_filter == "all" or row.get("tool") == args.tool_filter:
                agg_data[row["name"]] = row

    print(f"[{args.target}] {len(agg_data)} designs dans le CSV agrégé", flush=True)

    # BindCraft : PDB plats dans Accepted/ (pas de sous-dossiers ColabFold)
    pdb_map = {}
    for fname in os.listdir(args.colabfold_dir):
        if not fname.endswith(".pdb") or "Ranked" in fname:
            continue
        design = fname[:-4]
        pdb_map[design] = os.path.join(args.colabfold_dir, fname)

    print(f"[{args.target}] {len(pdb_map)} PDB rank_001 trouvés", flush=True)

    # Scorer
    results = []
    errors  = 0
    t_global = time.time()

    # Charger les designs déjà scorés pour reprise après timeout
    already_done = set()
    if os.path.exists(args.output):
        with open(args.output) as f_done:
            reader_done = csv.DictReader(f_done)
            for row_done in reader_done:
                already_done.add(row_done["name"])
        print(f"[{args.target}] {len(already_done)} designs déjà scorés — skip", flush=True)

    for idx, design_name in enumerate(sorted(agg_data.keys())):
        if design_name not in pdb_map:
            continue
        if design_name in already_done:
            continue
        pdb_path = pdb_map[design_name]
        t0 = time.time()
        try:
            scores = score_af2_complex(
                pdb_path, args.binder_chain, args.target_chain, sfxn, fr)

            # Fusionner avec données agrégées AF2
            agg = agg_data.get(design_name, {})
            row = {
                "name":    design_name,
                "target":  args.target,
                "tool":    agg.get("tool", "rfdiff_mpnn"),
                # Métriques AF2 depuis CSV agrégé
                "iptm":    agg.get("iptm", ""),
                "ipsae":   agg.get("ipsae", ""),
                "i_pae":   agg.get("i_pae", ""),
                "ptm":     agg.get("ptm", ""),
                "mean_plddt": agg.get("mean_plddt", ""),
                # Métriques PyRosetta
                **scores,
            }
            results.append(row)

            # Écriture incrémentale — permet reprise après timeout
            write_header = not os.path.exists(args.output) or idx == 0
            with open(args.output, "a" if not write_header else "w", newline="") as f_out:
                writer = csv.DictWriter(f_out, fieldnames=row.keys())
                if write_header:
                    writer.writeheader()
                writer.writerow(row)

            elapsed = time.time() - t0
            iptm_str = f"iPTM={float(agg['iptm']):.3f}" if agg.get("iptm") else "iPTM=?"
            print(f"[{idx+1:3d}/{len(agg_data)}] {design_name:40s} "
                  f"dG={scores['dG']:7.1f} SC={scores['SC']:.3f} "
                  f"dSASA={scores['dSASA']:6.0f} {iptm_str} {elapsed:.1f}s",
                  flush=True)

        except Exception as e:
            print(f"  ⚠ Erreur {design_name}: {e}", flush=True)
            errors += 1

    # Sauvegarde déjà faite ligne par ligne durant le run

    total = time.time() - t_global
    n = max(len(results), 1)
    print(f"\n=== RÉSUMÉ [{args.target}] ===", flush=True)
    print(f"Scorés  : {len(results)}", flush=True)
    print(f"Erreurs : {errors}", flush=True)
    print(f"Temps   : {total:.0f}s ({total/n:.1f}s/design)", flush=True)
    print(f"Output  : {args.output}", flush=True)

    # Corrélation PyRosetta vs iPTM
    paired = [(float(r["iptm"]), r["dG"], r["SC"], r["dSASA"])
              for r in results if r["iptm"]]
    if len(paired) >= 5:
        def pearson(x, y):
            n = len(x)
            mx, my = sum(x)/n, sum(y)/n
            num = sum((xi-mx)*(yi-my) for xi,yi in zip(x,y))
            den = (sum((xi-mx)**2 for xi in x)*sum((yi-my)**2 for yi in y))**0.5
            return num/den if den else 0

        iptms = [p[0] for p in paired]
        dgs   = [p[1] for p in paired]
        scs   = [p[2] for p in paired]
        dsasas= [p[3] for p in paired]

        print(f"\n=== CORRÉLATION PyRosetta(AF2) vs iPTM ({len(paired)} designs) ===",
              flush=True)
        print(f"r(dG,   iPTM) = {pearson(dgs,    iptms):.3f}  (attendu négatif)",
              flush=True)
        print(f"r(SC,   iPTM) = {pearson(scs,    iptms):.3f}  (attendu positif)",
              flush=True)
        print(f"r(dSASA,iPTM) = {pearson(dsasas, iptms):.3f}  (attendu positif)",
              flush=True)

        # Simulation filtrage
        threshold = 0.85
        n_good = sum(1 for p in paired if p[0] >= threshold)
        # Filtre composite : dG < -15 ET SC > 0.55
        filtered = [p for p in paired if p[1] < -15 and p[2] > 0.55]
        n_good_filtered = sum(1 for p in filtered if p[0] >= threshold)

        print(f"\n=== SIMULATION FILTRAGE (dG<-15 ET SC>0.55) ===", flush=True)
        print(f"Total               : {len(paired)}", flush=True)
        print(f"Bons (iPTM>{threshold}) : {n_good} ({100*n_good/len(paired):.0f}%)",
              flush=True)
        print(f"Après filtre        : {len(filtered)} designs", flush=True)
        if filtered:
            recall = n_good_filtered/max(n_good,1)
            precision = n_good_filtered/max(len(filtered),1)
            print(f"Bons gardés         : {n_good_filtered}", flush=True)
            print(f"Recall              : {recall:.2f} ({recall*100:.0f}% des bons conservés)",
                  flush=True)
            print(f"Precision           : {precision:.2f} ({precision*100:.0f}% des filtrés sont bons)",
                  flush=True)

if __name__ == "__main__":
    main()
