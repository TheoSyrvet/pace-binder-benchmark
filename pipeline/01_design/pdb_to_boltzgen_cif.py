#!/usr/bin/env python3
"""PDB croppé -> CIF BoltzGen avec _entity_poly_seq. Usage: script in.pdb out.cif"""
import sys
import gemmi

pdb_path, cif_path = sys.argv[1], sys.argv[2]

st = gemmi.read_structure(pdb_path)
st.setup_entities()
st.assign_label_seq_id()

# forcer le type polymere sur les entities (sinon pas de _entity_poly_seq)
for ent in st.entities:
    if ent.entity_type == gemmi.EntityType.Unknown:
        ent.entity_type = gemmi.EntityType.Polymer
    if ent.polymer_type == gemmi.PolymerType.Unknown:
        ent.polymer_type = gemmi.PolymerType.PeptideL

for chain in st[0]:
    seqids = [r.seqid.num for r in chain]
    print(f"  chaine {chain.name}: {len(chain)} res, range {min(seqids)}-{max(seqids)}")

doc = st.make_mmcif_document()
doc.write_file(cif_path)
print(f"  -> ecrit {cif_path}")

block = gemmi.cif.read(cif_path).sole_block()
eps = block.find_loop("_entity_poly_seq.num")
n = len(list(eps)) if eps else 0
print(f"  _entity_poly_seq : {'PRESENT (' + str(n) + ' lignes)' if n else 'ABSENT'}")
