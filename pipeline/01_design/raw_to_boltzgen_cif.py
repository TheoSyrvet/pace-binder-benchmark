#!/usr/bin/env python3
"""CIF RCSB raw -> CIF BoltzGen (chaine+segment), conserve _entity_poly_seq.
Usage: script raw.cif CHAIN START END out.cif"""
import sys, gemmi

raw, chain, start, end, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]

st = gemmi.read_structure(raw)
st.setup_entities()

sel = gemmi.Selection(f"/1/{chain}/{start}-{end}")
st2 = sel.copy_structure_selection(st)
st2.setup_entities()
st2.assign_label_seq_id()

for ch in st2[0]:
    ids = [r.seqid.num for r in ch]
    print(f"  chaine {ch.name}: {len(ch)} res, range {min(ids)}-{max(ids)}")

doc = st2.make_mmcif_document()
doc.write_file(out)
print(f"  -> {out}")

blk = gemmi.cif.read(out).sole_block()
eps = blk.find_loop("_entity_poly_seq.num")
n = len(list(eps)) if eps else 0
print(f"  _entity_poly_seq : {'PRESENT ('+str(n)+')' if n else 'ABSENT'}")
