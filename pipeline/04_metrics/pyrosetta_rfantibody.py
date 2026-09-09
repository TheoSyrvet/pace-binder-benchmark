#!/usr/bin/env python3
"""
PyRosetta scorer RFantibody (niveau 2, parallele OF3).
Tourne sur les memes 100 designs/cible qu'OF3 (pdbs/<cible>/<tag>.pdb, structures RF2).
Scoring generique dG/SC/dSASA, interface H_T (H=binder, T=target).
Fusionne metriques RF2 du top100 tsv. Pas d'ipSAE. Reprise incrementale.
"""
import os, sys, csv, time, argparse
import pyrosetta
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
from pyrosetta.rosetta.core.pose import Pose

BASE   = os.path.expanduser("~/master_thesis/results/rfantibody/rfantibody_top100")
PDBDIR = f"{BASE}/pdbs"
TSV    = f"{BASE}/ALL_top100.tsv"
OUTDIR = os.path.expanduser("~/master_thesis/results/ranking/rfantibody_pyrosetta")
os.makedirs(OUTDIR, exist_ok=True)
BINDER_CHAIN = "H"
TARGET_CHAIN = "T"

import tempfile

def clean_collapsed_residues(pdb_path):
    """Retire les residus dont tous les atomes ont les memes coords (placeholders
    RFantibody non resolus, surtout dans le target). Ecrit un PDB temporaire propre.
    Retourne (chemin_temp, n_retires)."""
    from collections import defaultdict, OrderedDict
    res_coords = OrderedDict()
    for l in open(pdb_path):
        if l.startswith("ATOM"):
            key = (l[21], l[22:26].strip())
            res_coords.setdefault(key, []).append((l[30:38], l[38:46], l[46:54]))
    collapsed = {k for k, v in res_coords.items() if len(set(v)) == 1 and len(v) > 1}
    if not collapsed:
        return pdb_path, 0
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False)
    for l in open(pdb_path):
        if l.startswith(("ATOM", "HETATM")):
            key = (l[21], l[22:26].strip())
            if key in collapsed:
                continue
        tmp.write(l)
    tmp.close()
    return tmp.name, len(collapsed)


def score_complex(pdb_path, binder_chain, target_chain, sfxn, fr):
    clean_path, n_removed = clean_collapsed_residues(pdb_path)
    pose = pose_from_pdb(clean_path)
    if clean_path != pdb_path:
        os.remove(clean_path)
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
    pyrosetta.init("-mute all -ex1 -ex2aro")
    sfxn = get_score_function(True)
    fr = FastRelax(sfxn, 1)
    rf2 = {}
    with open(TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["target"] == tgt:
                rf2[row["tag"]] = row
    pdir = f"{PDBDIR}/{tgt}"
    if not os.path.isdir(pdir):
        print(f"[{tgt}] dossier PDB absent — abort"); sys.exit(1)
    pdbs = {fn[:-4]: os.path.join(pdir, fn) for fn in os.listdir(pdir) if fn.endswith(".pdb")}
    print(f"[{tgt}] {len(pdbs)} PDB | {len(rf2)} lignes RF2", flush=True)
    out = f"{OUTDIR}/pyrosetta_{tgt}.csv"
    done = set()
    if os.path.exists(out):
        with open(out) as f:
            done = {r["name"] for r in csv.DictReader(f)}
        print(f"[{tgt}] {len(done)} deja scores — skip", flush=True)
    results, errors, t0 = [], 0, time.time()
    for idx, tag in enumerate(sorted(pdbs.keys())):
        if tag in done:
            continue
        try:
            sc = score_complex(pdbs[tag], BINDER_CHAIN, TARGET_CHAIN, sfxn, fr)
            r = rf2.get(tag, {})
            row = {
                "name": tag, "target": tgt, "tool": "rfantibody",
                "interaction_pae": r.get("interaction_pae", ""),
                "target_aligned_antibody_rmsd": r.get("target_aligned_antibody_rmsd", ""),
                "pred_lddt": r.get("pred_lddt", ""),
                "passe_gate": r.get("passe_gate", ""),
                **sc,
            }
            results.append(row)
            wh = not os.path.exists(out)
            with open(out, "a", newline="") as fo:
                w = csv.DictWriter(fo, fieldnames=row.keys())
                if wh: w.writeheader()
                w.writerow(row)
            print(f"[{idx+1:3d}/100] {tag:45s} dG={sc['dG']:7.1f} SC={sc['SC']:.3f} dSASA={sc['dSASA']:6.0f}", flush=True)
        except Exception as e:
            print(f"  WARN {tag}: {e}", flush=True); errors += 1
    dt = time.time() - t0
    print(f"\n=== [{tgt}] scores={len(results)} erreurs={errors} temps={dt:.0f}s ({dt/max(len(results),1):.1f}s/design) ===", flush=True)
    print(f"Output: {out}", flush=True)

if __name__ == "__main__":
    main()
