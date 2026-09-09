#!/usr/bin/env python3
"""Pour chaque cible : top-100 designs par iptm (rank_001) avec toutes les métriques,
+ moyenne des métriques par cible. Sort un CSV détail par cible + un CSV récap global."""
import json, glob, os, csv
import numpy as np

TARGETS = ["mdm2","pdl1","spike","tnfa","egfr","insulinr","vegf","stat3","kit","glp1r",
           "pdgfrb","fgfr2","trka","ha_head","ha_stem","il7ra","clathrin","cd28","ccr5","cxcr4"]
BINDER_LEN = 100
OUT = os.path.expanduser("~/master_thesis/results/ranking/ppdiff")
os.makedirs(OUT, exist_ok=True)

def metrics_for_design(d):
    js = glob.glob(f"{d}/*scores_rank_001*.json")
    if not js: return None
    try: s = json.load(open(js[0]))
    except: return None
    iptm, ptm = s.get("iptm"), s.get("ptm")
    plddt = np.array(s.get("plddt", [])); pae = np.array(s.get("pae", []))
    if plddt.size==0 or pae.size==0 or iptm is None: return None
    binder_plddt = float(plddt[:BINDER_LEN].mean())
    complex_plddt = float(plddt.mean())
    b=slice(0,BINDER_LEN); t=slice(BINDER_LEN,pae.shape[0])
    i_pae = float(np.concatenate([pae[b,t].ravel(), pae[t,b].ravel()]).mean())
    return dict(design=os.path.basename(d), iptm=round(iptm,4), ptm=round(ptm,4),
                binder_plddt=round(binder_plddt,3), complex_plddt=round(complex_plddt,3),
                i_pae=round(i_pae,3), max_pae=round(float(s.get("max_pae",0)),3))

recap = []
FIELDS = ["design","iptm","ptm","binder_plddt","complex_plddt","i_pae","max_pae"]
for tgt in TARGETS:
    base = os.path.expanduser(f"~/master_thesis/results/colabfold/ppdiff/{tgt}")
    rows = [m for d in sorted(glob.glob(f"{base}/{tgt}_ppdiff_*")) if (m:=metrics_for_design(d))]
    if not rows:
        print(f"{tgt}: aucun design lu"); continue
    rows.sort(key=lambda r: r["iptm"], reverse=True)
    top = rows[:100]
    # CSV détail par cible
    with open(f"{OUT}/{tgt}_top100.csv","w",newline="") as fh:
        w=csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader(); w.writerows(top)
    # moyennes
    means = {k: round(float(np.mean([r[k] for r in top])),4) for k in FIELDS if k!="design"}
    means["target"]=tgt; means["n_designs_total"]=len(rows); means["n_in_top"]=len(top)
    recap.append(means)
    print(f"{tgt:10} top{len(top):3} | iptm={means['iptm']:.3f} ptm={means['ptm']:.3f} "
          f"b_plddt={means['binder_plddt']:.2f} i_pae={means['i_pae']:.2f}")

# CSV récap global (moyennes par cible)
RFIELDS=["target","n_designs_total","n_in_top","iptm","ptm","binder_plddt","complex_plddt","i_pae","max_pae"]
with open(f"{OUT}/RECAP_means_per_target.csv","w",newline="") as fh:
    w=csv.DictWriter(fh, fieldnames=RFIELDS); w.writeheader()
    for r in recap: w.writerow({k:r.get(k) for k in RFIELDS})
print(f"\n>>> {OUT}/RECAP_means_per_target.csv ({len(recap)} cibles)")
print(f">>> {OUT}/<cible>_top100.csv (détail par cible)")
