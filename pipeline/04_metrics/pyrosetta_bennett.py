#!/usr/bin/env python3
"""
PyRosetta scorer Bennett (RFdiffusion+SolMPNN), niveau 2, parallele OF3.
Tourne sur les MEMES top100/cible qu'OF3 (par pae_interaction).
Structures : PDB af2_out/<design>.pdb. ATTENTION Bennett : chaine A=BINDER, B=TARGET.
Interface A_B. Reprise incrementale. 3 trimeres exclus.
"""
import os, sys, csv, time, glob, argparse
import pyrosetta
from pyrosetta import pose_from_file
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pose import Pose

TOP_DIR = os.path.expanduser("~/master_thesis/results/ranking/bennett_top100")
OUTDIR  = os.path.expanduser("~/master_thesis/results/ranking/bennett_pyrosetta")
os.makedirs(OUTDIR, exist_ok=True)
# Bennett : A=binder, B=target (inverse de boltzgen)
BINDER_CHAIN = "A"
TARGET_CHAIN = "B"

def score_complex(pdb_path, binder_chain, target_chain, sfxn, fr):
    pose = pose_from_file(pdb_path)
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

    topfile = f"{TOP_DIR}/{tgt}_top.txt"
    if not os.path.isfile(topfile):
        print(f"[{tgt}] pas de top.txt — abort"); sys.exit(1)
    designs = []  # (desc, pae, pdb)
    for line in open(topfile):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3: designs.append((p[0], p[1], p[2]))
    if not designs:
        print(f"[{tgt}] top.txt vide — abort"); sys.exit(1)

    pyrosetta.init("-mute all -ex1 -ex2aro")
    sfxn = get_score_function(True)
    fr = FastRelax(sfxn, 1)

    out = f"{OUTDIR}/pyrosetta_{tgt}.csv"
    done = set()
    if os.path.exists(out):
        with open(out) as f:
            done = {r["name"] for r in csv.DictReader(f)}
        print(f"[{tgt}] {len(done)} deja scores — skip", flush=True)

    print(f"[{tgt}] {len(designs)} designs top100", flush=True)
    results, errors, t0 = [], 0, time.time()
    for idx, (desc, pae, pdb) in enumerate(designs):
        if desc in done: continue
        if not os.path.isfile(pdb):
            print(f"  WARN {desc}: pdb absent", flush=True); errors += 1; continue
        try:
            sc = score_complex(pdb, BINDER_CHAIN, TARGET_CHAIN, sfxn, fr)
            row = {"name": desc, "target": tgt, "tool": "bennett",
                   "pae_interaction": pae, **sc}
            results.append(row)
            wh = not os.path.exists(out)
            with open(out, "a", newline="") as fo:
                w = csv.DictWriter(fo, fieldnames=row.keys())
                if wh: w.writeheader()
                w.writerow(row)
            print(f"[{idx+1:3d}/100] {desc:35s} dG={sc['dG']:7.1f} SC={sc['SC']:.3f} dSASA={sc['dSASA']:6.0f}", flush=True)
        except Exception as e:
            print(f"  WARN {desc}: {e}", flush=True); errors += 1
    dt = time.time() - t0
    print(f"\n=== [{tgt}] scores={len(results)} erreurs={errors} temps={dt:.0f}s ===", flush=True)

if __name__ == "__main__":
    main()
