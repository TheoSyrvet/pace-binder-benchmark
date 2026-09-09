#!/usr/bin/env python3
"""Ranking PPDiff : lit les scores ColabFold rank_001, calcule binder-plddt + i_pAE,
applique le filtre Step 1, sort les 100 meilleurs par iptm. Binder = chaine A (100 res)."""
import json, glob, os, sys
import numpy as np

TARGET = sys.argv[1]
BASE = os.path.expanduser(f"~/master_thesis/results/colabfold/ppdiff/{TARGET}")
BINDER_LEN = 100  # chaine A, longueur fixe PPDiff

# seuils Step 1 (validés)
TH = dict(iptm=0.5, ptm=0.55, binder_plddt=80.0, i_pae=10.85)

rows = []
design_dirs = sorted(glob.glob(f"{BASE}/{TARGET}_ppdiff_*"))
for d in design_dirs:
    js = glob.glob(f"{d}/*scores_rank_001*.json")
    if not js:
        continue
    try:
        s = json.load(open(js[0]))
    except Exception:
        continue
    iptm = s.get("iptm")
    ptm = s.get("ptm")
    plddt = np.array(s.get("plddt", []))
    pae = np.array(s.get("pae", []))
    if plddt.size == 0 or pae.size == 0 or iptm is None:
        continue
    # binder plddt = moyenne residus 0..BINDER_LEN-1
    binder_plddt = float(plddt[:BINDER_LEN].mean())
    # i_pAE = moyenne des blocs croises binder<->target
    b = slice(0, BINDER_LEN)
    t = slice(BINDER_LEN, pae.shape[0])
    cross = np.concatenate([pae[b, t].ravel(), pae[t, b].ravel()])
    i_pae = float(cross.mean())
    passed = (iptm > TH["iptm"] and ptm > TH["ptm"]
              and binder_plddt > TH["binder_plddt"] and i_pae < TH["i_pae"])
    rows.append(dict(design=os.path.basename(d), iptm=iptm, ptm=ptm,
                     binder_plddt=round(binder_plddt,2), i_pae=round(i_pae,2),
                     pass_step1=passed))

rows.sort(key=lambda r: r["iptm"], reverse=True)
n_pass = sum(r["pass_step1"] for r in rows)
print(f"=== {TARGET} : {len(rows)} designs lus, {n_pass} passent Step 1 ===")
print(f"{'design':<22} {'iptm':>5} {'ptm':>5} {'b_plddt':>8} {'i_pae':>6} {'pass':>5}")
for r in rows[:15]:
    print(f"{r['design']:<22} {r['iptm']:>5} {r['ptm']:>5} {r['binder_plddt']:>8} {r['i_pae']:>6} {str(r['pass_step1']):>5}")

# sauvegarde CSV complet + les 100 meilleurs (parmi ceux qui passent, sinon top iptm)
import csv
outdir = os.path.expanduser(f"~/master_thesis/results/ranking/ppdiff")
os.makedirs(outdir, exist_ok=True)
with open(f"{outdir}/{TARGET}_ranked.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["design","iptm","ptm","binder_plddt","i_pae","pass_step1"])
    w.writeheader(); w.writerows(rows)
passing = [r for r in rows if r["pass_step1"]]
top100 = (passing if len(passing)>=100 else rows)[:100]
with open(f"{outdir}/{TARGET}_top100.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["design","iptm","ptm","binder_plddt","i_pae","pass_step1"])
    w.writeheader(); w.writerows(top100)
print(f">>> sauvegardé: {outdir}/{TARGET}_ranked.csv ({len(rows)} designs) + {TARGET}_top100.csv ({len(top100)})")
