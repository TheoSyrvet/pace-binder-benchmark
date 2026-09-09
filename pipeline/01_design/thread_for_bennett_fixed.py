#!/usr/bin/env python3
"""Thread les sequences SolMPNN sur les backbones RFdiff pour Bennett.
Fix v3: predict.py attend CIBLE en chaine A (premiere) et BINDER en chaine B (seconde).
binderlen = splits[1].size() = taille binder = deuxieme chaine.
"""
import sys, os, re, glob
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain

aa1to3={'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN','E':'GLU','G':'GLY',
'H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE','P':'PRO','S':'SER','T':'THR',
'W':'TRP','Y':'TYR','V':'VAL','X':'GLY'}

target=sys.argv[1]
base=os.path.expanduser('~/master_thesis')
rfdir=f'{base}/results/rfdiffusion/final/{target}'
mpnndir=f'{base}/results/mpnn/final/{target}'
outdir=f'{base}/results/bennett/final/{target}/pdbs'
os.makedirs(outdir, exist_ok=True)

parser=PDBParser(QUIET=True)
n_written=0; n_skip=0

for fa in sorted(glob.glob(f'{mpnndir}/design_*/seqs/design_*.fa')):
    design=re.search(r'(design_\d+)', fa).group(1)
    bb=f'{rfdir}/{design}.pdb'
    if not os.path.exists(bb):
        n_skip+=1; continue
    lines=[l.rstrip('\n') for l in open(fa)]
    hdr=lines[0]
    native_seq=lines[1].strip()

    # Identifier le binder = chaine dont la seq native est all-G
    m=re.search(r"designed_chains=\[([^\]]*)\]", hdr)
    designed_chains=m.group(1).replace("'","").replace(" ","").split(',')
    segments=native_seq.split('/')
    binder_ch=None
    binder_seg_idx=None
    for idx,(ch,seg) in enumerate(zip(designed_chains, segments)):
        if len(seg)>0 and all(c=='G' for c in seg):
            binder_ch=ch
            binder_seg_idx=idx
            break
    if binder_ch is None:
        n_skip+=1; continue

    # Parser les samples
    samples=[]
    for i,l in enumerate(lines):
        if l.startswith('>') and 'sample=' in l:
            snum=re.search(r'sample=(\d+)', l).group(1)
            fullseq=lines[i+1].strip()
            parts=fullseq.split('/')
            if binder_seg_idx < len(parts):
                binder_seq=parts[binder_seg_idx].strip()
                samples.append((snum, binder_seq))

    model=parser.get_structure('x', bb)[0]
    chains=[c.id for c in model]
    target_chs=[c for c in chains if c!=binder_ch]

    for snum, binder_seq in samples:
        try:
            bchain=model[binder_ch]
        except KeyError:
            n_skip+=1; continue
        bres=list(bchain)
        if len(bres)!=len(binder_seq):
            n_skip+=1; continue
        for res,aa in zip(bres, binder_seq):
            res.resname=aa1to3.get(aa,'GLY')

        # predict.py attend: CIBLE->chaine A (premiere), BINDER->chaine B (seconde)
        # binderlen = splits[1].size() = taille binder
        new_s=Structure.Structure('x')
        new_m=Model.Model(0)
        new_s.add(new_m)

        # Cible en chaine A (premiere)
        new_tc=Chain.Chain('A')
        new_m.add(new_tc)
        rnum=1
        for tch in target_chs:
            if tch not in model:
                continue
            for res in model[tch]:
                rc=res.copy()
                rc.id=(rc.id[0], rnum, rc.id[2])
                new_tc.add(rc)
                rnum+=1

        # Binder en chaine B (seconde)
        new_bc=Chain.Chain('B')
        new_m.add(new_bc)
        for res in bres:
            new_bc.add(res.copy())

        outpdb=f'{outdir}/{target}_design_{design.split("_")[1]}_sample{snum}.pdb'
        io=PDBIO()
        io.set_structure(new_s)
        io.save(outpdb)
        n_written+=1

print(f'{target}: {n_written} PDB ecrits, {n_skip} skip')
