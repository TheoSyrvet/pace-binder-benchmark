#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess, csv, sys
from collections import defaultdict

INDEX_TO_TARGET = {
    1:"mdm2",2:"pdl1",3:"spike",4:"tnfa",5:"egfr",6:"insulinr",7:"vegf",
    8:"stat3",9:"kit",10:"glp1r",11:"pdgfrb",12:"fgfr2",13:"trka",
    14:"ha_head",15:"ha_stem",16:"il7ra",17:"clathrin",18:"cd28",19:"ccr5",20:"cxcr4",
}
TARGETS=[INDEX_TO_TARGET[i] for i in range(1,21)]
RFD_NORM={"hastem":"ha_stem","ha_stem":"ha_stem","ha_head":"ha_head","hahead":"ha_head"}

# (modele_affiche, composante, job_name, starttime, ventilation)
# RFdiffusion regroupe la generation RFdiff + le scoring Bennett (meme pipeline)
COMPONENTS=[
    ("RFdiffusion","rfdiff_gen","rfd-","2026-06-19","rfdname"),
    ("RFdiffusion","bennett_score","bennett_off","2026-06-19T17:00:00","array"),
    ("RFdiffusion","bennett_af2","bennett_af2","2026-06-25","array"),
    ("BindCraft","design","bc_final","2026-03-01","array"),
    ("BoltzGen","generation","boltzgen_final","2026-06-09","array"),
    ("BoltzGen","rerank","bg_rerank","2026-06-09","array"),
    ("PPDiff","gen_comp","ppdiff_comp","2026-06-11","array"),
    ("PPDiff","gen_final","ppdiff_final","2026-06-11","array"),
    ("PPDiff","scoring_cf","cf_ppd_clean","2026-06-11","array"),
    ("RFantibody","robust_genscore","rfab_robust","2026-06-12","array"),
    ("RFantibody","split_score","rfab_split","2026-06-12","array"),
    ("RFantibody","kit_mpnn","kit_mpnn","2026-06-12","array"),
]
VALID_STATES={"COMPLETED","TIMEOUT"}

def e2h(e):
    e=e.strip()
    if not e or e in ("INVALID","UNLIMITED"): return 0.0
    days=0
    if "-" in e: d,e=e.split("-",1); days=int(d)
    p=e.split(":")
    if len(p)==3: h,m,s=p
    elif len(p)==2: h,m,s="0",p[0],p[1]
    else: return 0.0
    try: return days*24+int(h)+int(m)/60.0+int(s)/3600.0
    except ValueError: return 0.0

def tgt_array(jid):
    b=jid.split(".")[0]
    if "_" not in b: return None
    i=b.split("_",1)[1]
    return INDEX_TO_TARGET.get(int(i)) if i.isdigit() else None

def tgt_rfd(jn):
    c=jn[len("rfd-"):]
    for suf in ("-final","-fix"):
        if c.endswith(suf): c=c[:-len(suf)]; break
    c=c.lower(); c=RFD_NORM.get(c,c)
    return c if c in TARGETS else None

def sacct(name,st):
    cmd=["sacct","-u","syrvet4","--starttime",st,"--name",name,
         "--format","JobID,JobName,Elapsed,State,Partition","--noheader","--parsable2"]
    try: return subprocess.run(cmd,capture_output=True,text=True,timeout=120).stdout.splitlines()
    except Exception as ex: print(f"[ERR {name}] {ex}",file=sys.stderr); return []

gpu=defaultdict(lambda:defaultdict(float)); jc=defaultdict(int); un=defaultdict(float)

def coll_rfd(st):
    cmd=["sacct","-u","syrvet4","--starttime",st,
         "--format","JobID,JobName,Elapsed,State,Partition","--noheader","--parsable2"]
    for line in subprocess.run(cmd,capture_output=True,text=True,timeout=180).stdout.splitlines():
        f=line.split("|")
        if len(f)<4: continue
        jid,jn,el,stt=f[0],f[1],f[2],f[3]
        if "." in jid or not jn.startswith("rfd-") or stt not in VALID_STATES: continue
        h=e2h(el); c=tgt_rfd(jn); jc[("RFdiffusion","rfdiff_gen")]+=1
        if c: gpu[("RFdiffusion","rfdiff_gen")][c]+=h
        else: un[("RFdiffusion","rfdiff_gen")]+=h

def coll_arr(mod,comp,name,st):
    for line in sacct(name,st):
        f=line.split("|")
        if len(f)<4: continue
        jid,jn,el,stt=f[0],f[1],f[2],f[3]
        part=f[4] if len(f)>4 else ""
        if "." in jid or jn!=name or stt not in VALID_STATES: continue
        if name=="bennett_off" and "cpu" in part.lower(): continue
        h=e2h(el); c=tgt_array(jid); jc[(mod,comp)]+=1
        if c: gpu[(mod,comp)][c]+=h
        else: un[(mod,comp)]+=h

print("=== Collecte GPU-h (COMPLETED+TIMEOUT) ===\n")
for mod,comp,name,st,vent in COMPONENTS:
    print(f"  {mod:12s}/{comp:16s} ({name}) >= {st}")
    if vent=="rfdname": coll_rfd(st)
    else: coll_arr(mod,comp,name,st)

with open("gpu_hours_long.csv","w",newline="") as fh:
    w=csv.writer(fh); w.writerow(["modele","composante","cible","gpu_h"])
    for (mod,comp),pt in sorted(gpu.items()):
        for c in TARGETS:
            if pt.get(c,0)>0: w.writerow([mod,comp,c,f"{pt[c]:.2f}"])

# Matrice par MODELE (somme des composantes par cible)
modeles=["RFdiffusion","BindCraft","BoltzGen","PPDiff","RFantibody"]
per_model_target=defaultdict(lambda:defaultdict(float))
for (mod,comp),pt in gpu.items():
    for c,h in pt.items(): per_model_target[mod][c]+=h

with open("gpu_hours_par_cible.csv","w",newline="") as fh:
    w=csv.writer(fh)
    w.writerow(["cible"]+modeles+["TOTAL"])
    coltot=defaultdict(float); grand=0.0
    for c in TARGETS:
        row=[c]; rt=0.0
        for mod in modeles:
            h=per_model_target[mod].get(c,0.0); row.append(f"{h:.1f}"); coltot[mod]+=h; rt+=h
        row.append(f"{rt:.1f}"); grand+=rt; w.writerow(row)
    w.writerow(["TOTAL"]+[f"{coltot[m]:.1f}" for m in modeles]+[f"{grand:.1f}"])

print("\n=== GPU-h par CIBLE et par MODELE ===")
hdr=f"{'cible':10s}"+"".join(f"{m[:11]:>12s}" for m in modeles)+f"{'TOTAL':>10s}"
print(hdr); print("-"*len(hdr))
coltot=defaultdict(float); grand=0.0
for c in TARGETS:
    line=f"{c:10s}"; rt=0.0
    for mod in modeles:
        h=per_model_target[mod].get(c,0.0); line+=f"{h:12.1f}"; coltot[mod]+=h; rt+=h
    line+=f"{rt:10.1f}"; grand+=rt; print(line)
print("-"*len(hdr))
tot=f"{'TOTAL':10s}"+"".join(f"{coltot[m]:12.1f}" for m in modeles)+f"{grand:10.1f}"
print(tot)

print("\n=== Detail composantes (verif) ===")
for (mod,comp),pt in sorted(gpu.items()):
    h=sum(pt.values()); n=jc.get((mod,comp),0); u=un.get((mod,comp),0.0)
    ex=f"  [!] {u:.1f}h NON attribuees" if u>0 else ""
    print(f"  {mod:12s}/{comp:16s}: {h:8.1f}h ({n} jobs){ex}")

print("\nFichiers: gpu_hours_par_cible.csv, gpu_hours_long.csv")
