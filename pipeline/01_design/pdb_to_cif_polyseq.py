#!/usr/bin/env python3
"""PDB croppe -> CIF BoltzGen avec _entity_poly_seq construit a la main.
Part du PDB (verite RFdiffusion), garde chaine+numerotation telles quelles.
Usage: script in.pdb out.cif"""
import sys, gemmi

pdb_path, cif_path = sys.argv[1], sys.argv[2]

st = gemmi.read_structure(pdb_path)
st.setup_entities()
st.assign_label_seq_id()

for ch in st[0]:
    ids = [r.seqid.num for r in ch]
    print(f"  chaine {ch.name}: {len(ch)} res, range {min(ids)}-{max(ids)}")

doc = st.make_mmcif_document()
block = doc.sole_block()

# Construire _entity_poly_seq manuellement depuis les residus
# (entity_id, num sequentiel, mon_id = nom du residu 3 lettres)
rows = []
for ent in st.entities:
    eid = ent.name
    for ch in st[0]:
        if ch.name in ent.subchains or ch.name == eid:
            for i, r in enumerate(ch, start=1):
                rows.append((eid, str(i), r.name))
            break

if rows:
    loop = block.init_loop("_entity_poly_seq.", ["entity_id", "num", "mon_id"])
    for row in rows:
        loop.add_row(list(row))

doc.write_file(cif_path)
print(f"  -> {cif_path}")

blk = gemmi.cif.read(cif_path).sole_block()
eps = blk.find_loop("_entity_poly_seq.num")
n = len(list(eps)) if eps else 0
print(f"  _entity_poly_seq : {'PRESENT ('+str(n)+')' if n else 'ABSENT'}")
