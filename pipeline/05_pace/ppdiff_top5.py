#!/usr/bin/env python3
import json, glob, os
import numpy as np
TARGETS = ["mdm2","pdl1","spike","tnfa","egfr","insulinr","vegf","stat3","kit","glp1r",
           "pdgfrb","fgfr2","trka","ha_head","ha_stem","il7ra","clathrin","cd28","ccr5","cxcr4"]
BINDER_LEN = 100
def m(d):
    js=glob.glob(f"{d}/*scores_rank_001*.json")
    if not js: return None
    try: s=json.load(open(js[0]))
    except: return None
    iptm,ptm=s.get("iptm"),s.get("ptm")
    pl=np.array(s.get("plddt",[])); pae=np.array(s.get("pae",[]))
    if pl.size==0 or pae.size==0 or iptm is None: return None
    bp=float(pl[:BINDER_LEN].mean())
    b=slice(0,BINDER_LEN); t=slice(BINDER_LEN,pae.shape[0])
    ip=float(np.concatenate([pae[b,t].ravel(),pae[t,b].ravel()]).mean())
    return (os.path.basename(d).split("_")[-1], round(iptm,3), round(ptm,3), round(bp,1), round(ip,2))
print(f"{'cible':<10}{'design':>7}{'iptm':>7}{'ptm':>7}{'b_plddt':>9}{'i_pae':>7}")
print("-"*47)
for tgt in TARGETS:
    base=os.path.expanduser(f"~/master_thesis/results/colabfold/ppdiff/{tgt}")
    rows=[r for d in sorted(glob.glob(f"{base}/{tgt}_ppdiff_*")) if (r:=m(d))]
    rows.sort(key=lambda x:x[1], reverse=True)
    for i,r in enumerate(rows[:5]):
        lbl=tgt if i==0 else ""
        print(f"{lbl:<10}{r[0]:>7}{r[1]:>7.3f}{r[2]:>7.3f}{r[3]:>9.1f}{r[4]:>7.2f}")
    print("-"*47)
