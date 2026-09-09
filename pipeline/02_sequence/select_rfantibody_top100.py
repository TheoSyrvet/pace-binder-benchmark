#!/usr/bin/env python3
"""
Sélection des 100 meilleurs designs RFantibody par cible, prêts pour OF3 -> PyRosetta.

Règle de sélection :
  1. Métrique-phare = interaction_pae (filtre RF2 officiel, Bennett et al. 2025, pAE<10).
  2. RMSD self-consistency = target_aligned_antibody_rmsd (col 4), critère "RF2-approved" (<2 A).
  3. Tri = rang combiné : rang(interaction_pae) + rang(target_aligned_rmsd), tous deux croissants.
  4. Gate pae<10 : si >=100 passent -> top 100 strict ; sinon gate relâché sur tout le pool.
"""

import os
import sys
import csv

BASE = os.path.expanduser("~/master_thesis/results/rfantibody/final")
OUTDIR = os.path.expanduser("~/master_thesis/results/rfantibody/rfantibody_top100")

TARGETS = ["mdm2", "pdl1", "spike", "tnfa", "egfr", "insulinr", "vegf", "stat3",
           "kit", "glp1r", "pdgfrb", "fgfr2", "trka", "ha_head", "ha_stem",
           "il7ra", "clathrin", "cd28", "ccr5", "cxcr4"]

PAE_GATE = 10.0
RMSD_COL = "target_aligned_antibody_rmsd"
PAE_COL  = "interaction_pae"
N_KEEP   = 100

EXPECTED_HEADER = [
    "interaction_pae", "pae", "pred_lddt",
    "target_aligned_antibody_rmsd", "target_aligned_cdr_rmsd",
    "framework_aligned_antibody_rmsd", "framework_aligned_cdr_rmsd",
    "framework_aligned_H1_rmsd", "framework_aligned_H2_rmsd", "framework_aligned_H3_rmsd",
    "framework_aligned_L1_rmsd", "framework_aligned_L2_rmsd", "framework_aligned_L3_rmsd",
    "tag",
]


def read_sc(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) != len(header):
                continue
            rows.append(dict(zip(header, parts)))
    return header, rows


def to_float(x):
    try:
        v = float(x)
        if v != v:
            return float("inf")
        return v
    except (ValueError, TypeError):
        return float("inf")


def select_target(target):
    sc = os.path.join(BASE, target, f"{target}_3_rf2.sc")
    if not os.path.isfile(sc):
        return [], {"target": target, "status": "SC_ABSENT", "n_total": 0,
                    "n_pass_gate": 0, "n_selected": 0, "gate_released": False}

    header, rows = read_sc(sc)

    if header != EXPECTED_HEADER:
        sys.stderr.write(f"[WARN] {target}: header inattendu, vérifie:\n"
                         f"  attendu: {EXPECTED_HEADER}\n  trouvé : {header}\n")

    for r in rows:
        r["_pae"]  = to_float(r.get(PAE_COL))
        r["_rmsd"] = to_float(r.get(RMSD_COL))
        r["_lddt"] = to_float(r.get("pred_lddt"))

    n_total = len(rows)
    passing = [r for r in rows if r["_pae"] < PAE_GATE]
    n_pass = len(passing)

    if n_pass >= N_KEEP:
        pool = passing
        gate_released = False
    else:
        pool = rows
        gate_released = True

    by_pae  = sorted(pool, key=lambda r: r["_pae"])
    rank_pae = {id(r): i for i, r in enumerate(by_pae)}
    by_rmsd = sorted(pool, key=lambda r: r["_rmsd"])
    rank_rmsd = {id(r): i for i, r in enumerate(by_rmsd)}

    for r in pool:
        r["_rank_pae"]  = rank_pae[id(r)]
        r["_rank_rmsd"] = rank_rmsd[id(r)]
        r["_rank_sum"]  = r["_rank_pae"] + r["_rank_rmsd"]
        r["_passe_gate"] = "1" if r["_pae"] < PAE_GATE else "0"

    pool_sorted = sorted(pool, key=lambda r: (r["_rank_sum"], r["_pae"]))
    selected = pool_sorted[:N_KEEP]

    info = {
        "target": target, "status": "OK", "n_total": n_total,
        "n_pass_gate": n_pass, "n_selected": len(selected),
        "gate_released": gate_released,
    }
    return selected, info


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    out_cols = ["target", "tag", "interaction_pae", "target_aligned_antibody_rmsd",
                "pred_lddt", "rank_pae", "rank_rmsd", "rank_sum", "passe_gate"]

    all_rows = []
    summary = []

    for t in TARGETS:
        selected, info = select_target(t)
        summary.append(info)

        per_target = os.path.join(OUTDIR, f"{t}_top100.tsv")
        with open(per_target, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(out_cols)
            for r in selected:
                row = [t, r.get("tag", ""),
                       f"{r['_pae']:.3f}", f"{r['_rmsd']:.3f}", f"{r['_lddt']:.3f}",
                       r["_rank_pae"], r["_rank_rmsd"], r["_rank_sum"], r["_passe_gate"]]
                w.writerow(row)
                all_rows.append(row)

        with open(os.path.join(OUTDIR, f"{t}_top100_tags.txt"), "w") as f:
            for r in selected:
                f.write(r.get("tag", "") + "\n")

    with open(os.path.join(OUTDIR, "ALL_top100.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(out_cols)
        w.writerows(all_rows)

    print(f"{'cible':<10} {'total':>6} {'pass_pae<10':>12} {'selected':>9} {'gate':>10}")
    print("-" * 52)
    tot_sel = 0
    for s in summary:
        gate = "RELACHE" if s["gate_released"] else "strict"
        print(f"{s['target']:<10} {s['n_total']:>6} {s['n_pass_gate']:>12} "
              f"{s['n_selected']:>9} {gate:>10}")
        tot_sel += s["n_selected"]
    print("-" * 52)
    print(f"TOTAL designs sélectionnés : {tot_sel}")
    print(f"\nSorties dans : {OUTDIR}")


if __name__ == "__main__":
    main()
