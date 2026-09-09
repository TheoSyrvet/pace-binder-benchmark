#!/usr/bin/env python3
"""Convertit les hotspots auth -> position sequentielle dans les YAML BoltzGen.
Lit le CIF, mappe chaque hotspot auth vers son rang dans la chaine cible.
Usage: script <cif> <chain> <hotspots_auth_csv>  -> imprime hotspots position"""
import sys, gemmi

cif, chain_id, hot_csv = sys.argv[1], sys.argv[2], sys.argv[3]
hotspots = [int(x) for x in hot_csv.split(",")]

st = gemmi.read_structure(cif)
model = st[0]

# trouver la chaine cible
chain = None
for ch in model:
    if ch.name == chain_id:
        chain = ch
        break
if chain is None:
    print(f"ERREUR chaine {chain_id} absente", file=sys.stderr)
    sys.exit(1)

# mapping auth -> position (rang 1-based dans la chaine)
auth_to_pos = {}
for i, res in enumerate(chain, start=1):
    auth_to_pos[res.seqid.num] = i

converted = []
for h in hotspots:
    if h in auth_to_pos:
        converted.append(auth_to_pos[h])
    else:
        print(f"  ATTENTION hotspot {h} absent de la chaine {chain_id}", file=sys.stderr)
        converted.append(None)

ok = [str(c) for c in converted if c is not None]
print(",".join(ok))
