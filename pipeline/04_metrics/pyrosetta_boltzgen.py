#!/usr/bin/env python3
"""
PyRosetta scorer BoltzGen (niveau 2, parallele OF3).
Tourne sur les MEMES top100/cible qu'OF3 (par final_rank natif).
Structures : .cif complexe dans intermediate_designs/<id>.cif (A=target, B=binder).
PyRosetta lit le cif directement (pose_from_file). Pas de residus collapsed.
Scoring generique dG/SC/dSASA, interface B_A. Fusionne metriques BoltzGen du CSV.
Reprise incrementale. kit exclu.
"""
import os, sys, csv, time, glob, argparse
import pyrosetta
from pyrosetta import pose_from_file
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
from pyrosetta.rosetta.core.pose import Pose

BG_DIR = os.path.expanduser("~/master_thesis/results/boltzgen/final")
OF3_IN = os.path.expanduser("~/master_thesis/shared/OpenFold3/inputs/boltzgen")
OUTDIR = os.path.expanduser("~/master_thesis/results/ranking/boltzgen_pyrosetta")
os.makedirs(OUTDIR, exist_ok=True)
BINDER_CHAIN = "B"
TARGET_CHAIN = "A"

def score_complex(cif_path, binder_chain, target_chain, sfxn, fr):
    pose = pose_from_file(cif_path)
    chain_map = {}
    for i in range(1, pose.num_chains() + 1):
        cid = pose.pdb_info().chain(pose.chain_begin(i))
        chain_map[cid] = i
    if binder_chain not in chain_map or target_chain not in chain_map:
        raise ValueError(f"chaines {binder_chain}/{target_chain} absentes (dispo: {list(chain_map.keys())})")
    tf = TaskFactory(); tf.push_back(RestrictToRepacking())
    PackRotamersMover(sfxn, tf.create_task_and_apply_taskoperations(pose)).apply(pose)
    fr.apply(pose)
    bnum = chain_map[binder_chain]
    bstart, bend = pose.chain_begin(bnum), pose.chain_end(bnum)
    blen = bend - bstart + 1
    binder_pose = Pose(pose, bstart, bend)
    binder_energy = sfxn(binder_pose)
    hydro = set("VILMFYW")
    bseq = "".join(pose.residue(i).name1() for i in range(bstart, bend + 1))
    surf_hydro = sum(1 for aa in bseq if aa in hydro) / max(len(bseq), 1)
    ia = InterfaceAnalyzerMover()
    ia.set_interface(f"{binder_chain}_{target_chain}")
    ia.set_scorefunction(sfxn)
    ia.set_compute_interface_energy(True)
    ia.set_compute_separated_sasa(True)
    ia.set_compute_interface_sc(True)
    ia.set_compute_interface_delta_hbond_unsat(True)
    ia.set_pack_input(True)
    ia.set_pack_separated(True)
    ia.apply(pose)
    dG = ia.get_interface_dG(); dSASA = ia.get_interface_delta_sasa()
    nres = ia.get_num_interface_residues(); unsat = ia.get_interface_delta_hbond_unsat()
    pstat = ia.get_interface_packstat(); hbe = ia.get_total_Hbond_E()
    cmplx = ia.get_complex_energy()
    scf = ShapeComplementarityFilter(); scf.jump_id(1)
    sc_val = scf.score(pose)
    return {
        "dG": round(dG, 3), "dSASA": round(dSASA, 1), "SC": round(sc_val, 4),
        "n_InterfaceResidues": nres, "n_InterfaceHbonds": round(hbe, 3),
        "n_InterfaceUnsatHbonds": round(unsat, 1), "PackStat": round(pstat, 4),
        "Binder_Energy_Score": round(binder_energy, 3),
        "Binder_Energy_per_res": round(binder_energy / max(blen, 1), 3),
        "Surface_Hydrophobicity": round(surf_hydro, 4),
        "dG_dSASA_ratio": round(dG / max(dSASA, 1), 6),
        "complex_energy": round(cmplx, 3),
        "binder_len": blen, "binder_seq": bseq,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    args = ap.parse_args()
    tgt = args.target

    # metriques BoltzGen natives (final_rank, design_to_target_iptm...) du CSV
    bg_csv = f"{BG_DIR}/{tgt}/final_ranked_designs/all_designs_metrics.csv"
    bgmet = {}
    if os.path.isfile(bg_csv):
        with open(bg_csv) as f:
            for r in csv.DictReader(f):
                bgmet[r["id"]] = r

    # les top100 = ceux qu'on a prepares pour OF3 (inputs/boltzgen/<tgt>/<tgt>__<id>.json)
    jsons = glob.glob(f"{OF3_IN}/{tgt}/{tgt}__*.json")
    ids = [os.path.basename(j)[:-5].split("__", 1)[1] for j in jsons]
    if not ids:
        print(f"[{tgt}] aucun design (pas de JSON OF3) — abort"); sys.exit(1)

    pyrosetta.init("-mute all -ex1 -ex2aro")
    sfxn = get_score_function(True)
    fr = FastRelax(sfxn, 1)

    out = f"{OUTDIR}/pyrosetta_{tgt}.csv"
    done = set()
    if os.path.exists(out):
        with open(out) as f:
            done = {r["name"] for r in csv.DictReader(f)}
        print(f"[{tgt}] {len(done)} deja scores — skip", flush=True)

    print(f"[{tgt}] {len(ids)} designs top100", flush=True)
    results, errors, t0 = [], 0, time.time()
    for idx, did in enumerate(sorted(ids)):
        if did in done:
            continue
        cif = f"{BG_DIR}/{tgt}/intermediate_designs/{did}.cif"
        if not os.path.isfile(cif):
            print(f"  WARN {did}: cif absent", flush=True); errors += 1; continue
        try:
            sc = score_complex(cif, BINDER_CHAIN, TARGET_CHAIN, sfxn, fr)
            m = bgmet.get(did, {})
            row = {
                "name": did, "target": tgt, "tool": "boltzgen",
                "final_rank": m.get("final_rank", ""),
                "design_to_target_iptm": m.get("design_to_target_iptm", ""),
                "min_design_to_target_pae": m.get("min_design_to_target_pae", ""),
                **sc,
            }
            results.append(row)
            wh = not os.path.exists(out)
            with open(out, "a", newline="") as fo:
                w = csv.DictWriter(fo, fieldnames=row.keys())
                if wh: w.writeheader()
                w.writerow(row)
            print(f"[{idx+1:3d}/100] {did:20s} dG={sc['dG']:7.1f} SC={sc['SC']:.3f} dSASA={sc['dSASA']:6.0f}", flush=True)
        except Exception as e:
            print(f"  WARN {did}: {e}", flush=True); errors += 1
    dt = time.time() - t0
    print(f"\n=== [{tgt}] scores={len(results)} erreurs={errors} temps={dt:.0f}s ({dt/max(len(results),1):.1f}s/design) ===", flush=True)
    print(f"Output: {out}", flush=True)

if __name__ == "__main__":
    main()
