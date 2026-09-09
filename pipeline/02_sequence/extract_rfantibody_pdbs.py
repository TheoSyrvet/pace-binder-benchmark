#!/usr/bin/env python3
"""
Extrait les PDB des 100 meilleurs designs RFantibody par cible depuis <cible>_3_rf2_full.qv.
"""
import os
import sys

BASE = os.path.expanduser("~/master_thesis/results/rfantibody/final")
TOPDIR = os.path.expanduser("~/master_thesis/results/rfantibody/rfantibody_top100")
PDBDIR = os.path.join(TOPDIR, "pdbs")

TARGETS = ["mdm2", "pdl1", "spike", "tnfa", "egfr", "insulinr", "vegf", "stat3",
           "kit", "glp1r", "pdgfrb", "fgfr2", "trka", "ha_head", "ha_stem",
           "il7ra", "clathrin", "cd28", "ccr5", "cxcr4"]


def load_wanted_tags(target):
    p = os.path.join(TOPDIR, f"{target}_top100_tags.txt")
    with open(p) as f:
        return set(line.strip() for line in f if line.strip())


def extract_from_qv(qv_path, wanted_tags):
    extracted = {}
    cur_tag = None
    cur_lines = None
    capturing = False
    with open(qv_path) as f:
        for line in f:
            if line.startswith("QV_TAG"):
                if capturing and cur_tag is not None:
                    extracted[cur_tag] = cur_lines
                parts = line.split()
                cur_tag = parts[1] if len(parts) > 1 else None
                if cur_tag in wanted_tags and cur_tag not in extracted:
                    capturing = True
                    cur_lines = []
                else:
                    capturing = False
                    cur_lines = None
            elif line.startswith("QV_SCORE"):
                continue
            else:
                if capturing:
                    cur_lines.append(line)
        if capturing and cur_tag is not None:
            extracted[cur_tag] = cur_lines
    return extracted


def main():
    os.makedirs(PDBDIR, exist_ok=True)
    grand_total = 0
    grand_missing = 0
    print(f"{'cible':<10} {'voulus':>7} {'extraits':>9} {'manquants':>10}")
    print("-" * 40)
    for t in TARGETS:
        qv = os.path.join(BASE, t, f"{t}_3_rf2_full.qv")
        if not os.path.isfile(qv):
            print(f"{t:<10} {'?':>7} {'QV_ABSENT':>9}")
            continue
        wanted = load_wanted_tags(t)
        outdir = os.path.join(PDBDIR, t)
        os.makedirs(outdir, exist_ok=True)
        extracted = extract_from_qv(qv, wanted)
        n_written = 0
        for tag, lines in extracted.items():
            outpath = os.path.join(outdir, f"{tag}.pdb")
            with open(outpath, "w") as f:
                f.writelines(lines)
            n_written += 1
        missing = wanted - set(extracted.keys())
        if missing:
            sys.stderr.write(f"[WARN] {t}: {len(missing)} tags introuvables:\n")
            for m in sorted(missing):
                sys.stderr.write(f"        {m}\n")
        grand_total += n_written
        grand_missing += len(missing)
        print(f"{t:<10} {len(wanted):>7} {n_written:>9} {len(missing):>10}")
    print("-" * 40)
    print(f"TOTAL extraits : {grand_total} / attendus {len(TARGETS)*100}, manquants : {grand_missing}")
    print(f"\nPDB dans : {PDBDIR}/<cible>/<tag>.pdb")


if __name__ == "__main__":
    main()
