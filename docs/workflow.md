# Protein Binder Design Pipeline — WORKFLOW

**Thesis:** Integrating RFdiffusion and Emerging Binder-Design Tools into an HPC Workflow
**Author:** Théo Syrvet — IGG, University of Geneva
**Supervisors:** Nicolas Hulo, Jose De Abreu Nunes
**Cluster:** Baobab HPC (SLURM) — account `falcone_teach_generic`, partition `shared-gpu` (GPU) / `shared-cpu` (CPU)
**Base directory:** `~/master_thesis/` (mirrored to `/acanas/biosc/master_thesis/`)

> This document supersedes the earlier WORKFLOW.md. The benchmark now covers **five** tools across **twenty** therapeutic targets, with a standardized ColabFold → PyRosetta evaluation chain.

---

## Overview

Five independent binder-design tools are benchmarked on each target. Each tool runs with its own optimal protocol. The standardized evaluation chain (AF2 confidence filter → PyRosetta energetic scoring) allows uniform comparison across tools.

```
Target PDB (cropped)
    │
    ├── RFdiffusion ──→ SolMPNN ─────────→ ColabFold ─┐
    ├── BindCraft ─────(internal AF2)────────────────┤
    ├── BoltzGen ──────(internal Boltz refold)───────┤
    ├── PPDiff ────────→ ColabFold ──────────────────┤
    └── RFantibody ────→ ColabFold ──────────────────┤
                                                       │
                                  AF2 filter (iPTM>0.5, pTM>0.55,
                                  pLDDT>0.8, i_pAE<10.85 Å)
                                                       │
                                  PyRosetta (dG<-15, SC>0.55, dSASA>800)
                                                       │
                                  per-target comparison CSV
```

**Target volume: ~800 designs per tool per target (100 backbones × 8 sequences for RFdiffusion; comparable for others).**

---

## Targets (20)

All cropped structures in `~/master_thesis/shared/targets/cropped_final/` (PDB) and `cropped_final_cif/` (CIF for BoltzGen).

| Target | PDB | Chain(s) | Segment | Hotspots |
|--------|-----|----------|---------|----------|
| mdm2 | 1YCR | A | 25–109 | 54,57,61,93,96,100 |
| pdl1 | 4ZQK | A | 19–127 | 54,56,113,115,121,123,124,125 |
| spike | 6M0J | E | 333–527 | 417,449,486,489,493,500,501 |
| tnfa | 2AZ5 | A+B | 10–157 | 57,59,119,120,151 |
| egfr | 1MOX | A | 1–190 | 45,50,88,89,101 |
| insulinr | 4ZXB | A | 470–910 | 479,484,498,704,708,710 |
| vegf | 1VPF | A+B | 14–107 | 17,21,25,79 |
| stat3 | 1BG1 | A | 580–680 | 591,609,611,612,613,636 |
| kit | 2E9W | A | 33–507 | 386,418,440,441,495,505 |
| glp1r | 5VAI | R | 29–145 | 39,67,69,99,101,121 |
| pdgfrb | 3MJG | A | 160–314 | 247,251,263,273 |
| fgfr2 | 1DJS | A | 140–284 | 207,241,251,257,281 |
| trka | 2IFG | A | 282–382 | 314,327,353 |
| ha_head | 4HMG | A | 50–260 | 98,136,153,183,190,194,195 |
| ha_stem | 3R2X | A+B | A11–324 / B1–171 | A18,38,40 / B21,41,45,49,52 |
| il7ra | 3DI3 | B | 17–202 | 119,124,152,191,192,194 |
| clathrin | 5M61 | A | 4–363 | 89,91,96,98 |
| cd28 | 1YJD | C | 1–118 | 99,100,101,103,105 |
| ccr5 | 4MBS | A | 19–313 | 168,172,177,180,190,262,264 |
| cxcr4 | 3ODU | A | 155–275 | 171,183,187,203,255,262 |

Same hotspots are imposed on every tool for a given target, ensuring binders aim at the same site (fair comparison).

---

## Containers & environments

| Tool | Mode | Location |
|------|------|----------|
| RFdiffusion | Apptainer | `RFdiffusion/containers/rfdiffusion_fixed.sif` (12 GB) |
| ProteinMPNN / SolMPNN | Apptainer | `shared/ProteinMPNN/container/pytorch.sif` (3 GB) |
| PPDiff | Apptainer | `PPDiff/ppdiff.sif` (5.3 GB) |
| RFantibody | Apptainer | `RFantibody/rfantibody.sif` (7.5 GB) |
| BindCraft | conda | env `BindCraft` (JAX + PyRosetta + AF2) |
| BoltzGen | conda | env `boltzgen` |
| ColabFold | module | `GCC/11.3.0 OpenMPI/4.1.4 ColabFold/1.5.2-CUDA-11.7.0` |

All Apptainer runs use `apptainer exec --nv` for GPU access.

---

## Tool 1 — RFdiffusion + SolMPNN  (reference pipeline)

**Strategy:** 100 backbones/target → 8 sequences/backbone = 800 sequences → ColabFold → AF2 filter → PyRosetta.

**RFdiffusion (per target):**
```bash
srun apptainer exec --nv \
  ~/master_thesis/RFdiffusion/containers/rfdiffusion_fixed.sif \
  python3.9 /app/RFdiffusion/scripts/run_inference.py \
  inference.input_pdb=~/master_thesis/shared/targets/cropped_final/1YCR_mdm2_final.pdb \
  'contigmap.contigs=[A25-109/0 50-150]' \
  'ppi.hotspot_res=[A54,A57,A61,A93,A96,A100]' \
  inference.output_prefix=~/master_thesis/results/rfdiffusion/final/mdm2/design \
  inference.model_directory_path=~/master_thesis/RFdiffusion/models \
  inference.schedule_directory_path=~/master_thesis/RFdiffusion/schedules \
  inference.num_designs=100 \
  denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0
```
Key points: `noise_scale=0` (deterministic, optimal for binders); hotspots required; IGSO3 schedule cached externally (`schedule_directory_path`).

**SolMPNN (sequence design):**
```bash
sbatch ~/master_thesis/shared/ProteinMPNN/scripts/run_solmpnn_final_array.sh
```
Uses soluble weights, `v_48_020`, `--num_seq_per_target 8 --sampling_temp 0.1`. Output: `results/mpnn/final/<target>/design_<N>/seqs/design_<N>.fa`.

---

## Tool 2 — BindCraft

**Strategy:** 100 accepted designs/target after internal AF2 filtering.
```bash
sbatch ~/master_thesis/BindCraft/scripts/run_bindcraft_final_array.sh
```
4-stage multimer protocol; soluble MPNN; cysteines excluded. **Critical:** `save_design_animations: false` (ffmpeg crash workaround). Internal AF2 metrics carry a known positive bias → BindCraft designs are re-scored with the common PyRosetta chain for fair comparison. No ColabFold needed.

---

## Tool 3 — BoltzGen

**Strategy:** generate 10000 designs, keep best 1000 by internal score.
```bash
sbatch ~/master_thesis/BoltzGen/scripts/run_boltzgen_final_array.sh
```
Input: YAML + **CIF** target (binder length 80–120). `--num_designs 10000 --budget 1000 --reuse --num_workers 0` (DataLoader deadlock fix). Internal Boltz refold metrics — bias as with BindCraft.

> **CIF caveat:** target CIF must contain the `_entity_poly_seq` block, otherwise BoltzGen fails parsing. Regenerate via gemmi from the original RCSB CIF (do NOT convert from PDB).

---

## Tool 4 — PPDiff

**Strategy:** 3 sequential runs × 267 designs = 801 sequences/target (no in-run resume; splitting limits loss on 12 h timeout).
```bash
sbatch ~/master_thesis/PPDiff/scripts/run_ppdiff_final_array.sh
```
fairseq compiled at job start; repaired checkpoint `checkpoint_fixed2.pt`. Binder uses dummy coords `[0,0,i]` (documented limitation). Outputs only sequences → require ColabFold.

> **Output-name caveat:** PPDiff appends a `1` suffix to output dirs (`run11`, `run21`, `run31`). Merge must read the suffixed dirs.

---

## Tool 5 — RFantibody  (nanobodies)

**Strategy:** 200 backbones × 4 sequences = 800 designs/target on a nanobody framework.
```bash
sbatch ~/master_thesis/RFantibody/run_rfantibody_final_array.sh
```
Pipeline: RFdiffusion_Ab (CDR diffusion, `-l "H1:6-8,H2:6-7,H3:8-14"`) → antibody ProteinMPNN (`-n 4 -t 0.2`) → RF2 (`-r 10`) → `qvscorefile`. Add `--no-trajectory` for spike/ha_stem/clathrin/cxcr4 (rotation-matrix bug). RF2 gives no iPTM → designs go through ColabFold.

---

## ColabFold validation  (RFdiffusion, PPDiff, RFantibody)

```bash
module load GCC/11.3.0 OpenMPI/4.1.4 ColabFold/1.5.2-CUDA-11.7.0
colabfold_batch --model-type alphafold2_multimer_v3 \
                --num-models 5 --num-recycle 3 <input.fasta> <out>
```
FASTA format: `>name\nBINDER_SEQ:TARGET_SEQ`. Use `alphafold2_multimer_v3` for everything (the `alphafold2_ptm` weights are corrupted on Baobab). A100 GPUs recommended.

Per-tool ColabFold launch scripts:
```bash
sbatch ~/master_thesis/shared/ColabFold/scripts/run_colabfold_final_array.sh        # RFdiffusion
sbatch ~/master_thesis/shared/ColabFold/scripts/run_colabfold_ppdiff_array.sh       # PPDiff
sbatch ~/master_thesis/shared/ColabFold/scripts/run_colabfold_rfantibody_array.sh   # RFantibody
```

---

## Evaluation chain (all tools)

**Step 1 — AF2 confidence filter** (mirrors BindCraft `default_filters.json`):
- iPTM > 0.5
- pTM > 0.55
- pLDDT > 0.8
- interface pAE < 10.85 Å  (= 0.35 × 31, BindCraft normalisation undone)

**Step 2 — PyRosetta energetic scoring** (CPU):
```bash
sbatch ~/master_thesis/shared/scripts/run_pyrosetta_array.sh
```
- dG < −15 REU
- shape complementarity > 0.55
- dSASA > 800 Å²

Script `pyrosetta_post_colabfold_v2.py` is incremental (writes line-by-line, skips already-scored designs → safe to resume after a 12 h timeout). BindCraft designs skip the AF2 filter (already internally filtered) but still go through PyRosetta.

---

## SLURM essentials

```bash
squeue -u $USER                       # queue
sbatch --array=3,4,17 script.sh       # relaunch selected targets only
scancel <jobid>
```
All partitions cap at **12 h**; every long step has a resume mechanism (skip already-done outputs). Always set `#SBATCH --account=falcone_teach_generic`.

---

## Quick test (RFdiffusion → SolMPNN → ColabFold, single target)

```bash
TARGET=mdm2
# 1. backbones
sbatch ~/master_thesis/RFdiffusion/scripts/final/run_rfdiff_${TARGET}_final.sh
# 2. sequences (after RFdiffusion done)
sbatch --array=1 ~/master_thesis/shared/ProteinMPNN/scripts/run_solmpnn_final_array.sh
# 3. fold & score
sbatch --array=1 ~/master_thesis/shared/ColabFold/scripts/run_colabfold_final_array.sh
# 4. PyRosetta
sbatch --array=1 ~/master_thesis/shared/scripts/run_pyrosetta_array.sh
```
