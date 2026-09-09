# Deployment notes

Per tool deployment recipes, parameters and known failure modes for the five design
tools in the benchmark. The evaluation method itself is described in the README.

Everything below assumes PACE_ROOT points at your working directory and that jobs
run on a SLURM cluster with GPU partitions. Replace the account name in the SLURM
directives with your own allocation.

## Targets

Twenty therapeutic targets, cropped to the domain carrying the binding site.
Cropped structures live under $PACE_ROOT/shared/targets/cropped_final/ as PDB and
under cropped_final_cif/ as CIF, the latter being required by BoltzGen.

The same hotspots are imposed on every tool for a given target, so that all methods
aim at the same geometric site. A poorly placed site would then penalise every
method equally on that target rather than favouring one of them.

| Target | PDB | Chains | Segment | Hotspots |
|--------|-----|--------|---------|----------|
| mdm2 | 1YCR | A | 25 to 109 | 54,57,61,93,96,100 |
| pdl1 | 4ZQK | A | 19 to 127 | 54,56,113,115,121,123,124,125 |
| spike | 6M0J | E | 333 to 527 | 417,449,486,489,493,500,501 |
| tnfa | 2AZ5 | A+B | 10 to 157 | 57,59,119,120,151 |
| egfr | 1MOX | A | 1 to 190 | 45,50,88,89,101 |
| insulinr | 4ZXB | A | 470 to 910 | 479,484,498,704,708,710 |
| vegf | 1VPF | A+B | 14 to 107 | 17,21,25,79 |
| stat3 | 1BG1 | A | 580 to 680 | 591,609,611,612,613,636 |
| kit | 2E9W | A | 33 to 507 | 386,418,440,441,495,505 |
| glp1r | 5VAI | R | 29 to 145 | 39,67,69,99,101,121 |
| pdgfrb | 3MJG | A | 160 to 314 | 247,251,263,273 |
| fgfr2 | 1DJS | A | 140 to 284 | 207,241,251,257,281 |
| trka | 2IFG | A | 282 to 382 | 314,327,353 |
| ha_head | 4HMG | A | 50 to 260 | 98,136,153,183,190,194,195 |
| ha_stem | 3R2X | A+B | A11 to 324, B1 to 171 | A18,38,40 and B21,41,45,49,52 |
| il7ra | 3DI3 | B | 17 to 202 | 119,124,152,191,192,194 |
| clathrin | 5M61 | A | 4 to 363 | 89,91,96,98 |
| cd28 | 1YJD | C | 1 to 118 | 99,100,101,103,105 |
| ccr5 | 4MBS | A | 19 to 313 | 168,172,177,180,190,262,264 |
| cxcr4 | 3ODU | A | 155 to 275 | 171,183,187,203,255,262 |

## Containers and environments

The five tools depend on mutually incompatible stacks, so each runs in its own
container or conda environment. Conda specifications are exported under env/conda.

| Tool | Mode | Location |
|------|------|----------|
| RFdiffusion | Apptainer | RFdiffusion/containers/rfdiffusion_fixed.sif, 12 GB |
| ProteinMPNN and SolMPNN | Apptainer | shared/ProteinMPNN/container/pytorch.sif, 3 GB |
| PPDiff | Apptainer | PPDiff/ppdiff.sif, 5.3 GB |
| RFantibody | Apptainer | RFantibody/rfantibody.sif, 7.5 GB |
| BindCraft | conda | env BindCraft, JAX with PyRosetta and AlphaFold2 |
| BoltzGen | conda | env boltzgen |
| OpenFold3 | conda | env openfold3 |
| MM-PBSA and trajectory analysis | conda | envs gmxmmpbsa and mdanalysis |
| ColabFold | module | GCC/11.3.0 OpenMPI/4.1.4 ColabFold/1.5.2-CUDA-11.7.0 |

Apptainer runs need apptainer exec --nv for GPU access.

## RFdiffusion with SolMPNN

One thousand backbones per target, each carrying a single sequence. That ratio is
not a preference but a constraint. When FastRelax is active, dl_binder_design
refuses to generate more than one sequence per backbone, and the guard is coded into
the pipeline rather than exposed as a setting. Keeping one relaxation cycle for the
quality it brings therefore means trading sequence diversity for structural
diversity, and reaching a comparable design count means generating more backbones.

Backbone generation:

```bash
srun apptainer exec --nv \
  $PACE_ROOT/RFdiffusion/containers/rfdiffusion_fixed.sif \
  python3.9 /app/RFdiffusion/scripts/run_inference.py \
  inference.input_pdb=$PACE_ROOT/shared/targets/cropped_final/1YCR_mdm2_final.pdb \
  'contigmap.contigs=[A25-109/0 50-150]' \
  'ppi.hotspot_res=[A54,A57,A61,A93,A96,A100]' \
  inference.output_prefix=$PACE_ROOT/results/rfdiffusion/final/mdm2/design \
  inference.model_directory_path=$PACE_ROOT/RFdiffusion/models \
  inference.schedule_directory_path=$PACE_ROOT/RFdiffusion/schedules \
  inference.num_designs=1000 \
  denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0
```

Noise scales are set to zero, which is deterministic and works better for binders.
Hotspots are required. The IGSO3 schedule is cached outside the container through
schedule_directory_path, otherwise it is recomputed at every launch.

Sequence design runs through dl_binder_design with soluble weights, model v_48_020
and a sampling temperature of 0.1, which keeps sampling conservative and close to
the model optimum. Note that the FASTA produced by SolMPNN contains the redesigned
binder only, not the fixed target sequence, which matters when the complex has to be
rebuilt downstream.

## BindCraft

One hundred accepted designs per target after internal AlphaFold2 filtering, using
the four stage multimer protocol with soluble MPNN and cysteines excluded.

Set save_design_animations: false. Leaving it on crashes the run through ffmpeg.

BindCraft filters its own designs with AlphaFold2 Multimer, the same model family it
optimises against, so its internal metrics carry a favourable bias. Its designs
therefore skip nothing and go through the common scoring chain like every other
tool. No ColabFold pass is needed since folding is already internal.

## BoltzGen

Generate a large pool and keep the best by internal score. Input is a YAML file plus
a CIF target, with binder length between 80 and 120 residues.

```
--num_designs 10000 --budget 1000 --reuse --num_workers 0
```

Two failure modes are worth knowing in advance.

The target CIF must contain the _entity_poly_seq block or BoltzGen fails while
parsing. Regenerate it with gemmi from the original RCSB CIF. Converting from PDB
does not produce that block.

Setting --num_workers 0 is not a performance choice. Any other value deadlocks the
data loader.

## PPDiff

Three sequential runs of 267 designs per target. PPDiff has no in-run resume, so
splitting the work limits what is lost when a job hits the twelve hour wall time.

fairseq is compiled at job start. The shipped checkpoint is corrupted and has to be
repaired, the working one being checkpoint_fixed2.pt. The binder is passed with
dummy coordinates [0,0,i], which is a documented limitation of the method rather
than a workaround.

PPDiff outputs sequences only, so its designs need a ColabFold pass before scoring.

One naming trap. PPDiff appends a 1 suffix to its output directories, producing
run11, run21 and run31 rather than run1, run2 and run3. Any merge step
has to read the suffixed names.

## RFantibody

Two hundred backbones with four sequences each, on a nanobody framework. The chain
is CDR diffusion, then antibody ProteinMPNN, then RoseTTAFold2 validation.

```
RFdiffusion_Ab   -l "H1:6-8,H2:6-7,H3:8-14"
ProteinMPNN_Ab   -n 4 -t 0.2
RF2              -r 10
```

Add --no-trajectory for spike, ha_stem, clathrin and cxcr4. Without it these four
targets hit a rotation matrix bug.

RoseTTAFold2 returns no ipTM, so RFantibody designs go through ColabFold like the
others before common scoring.

## ColabFold pass

Needed for RFdiffusion, PPDiff and RFantibody, whose native pipelines either produce
sequences only or use a predictor that does not report interface confidence.

```bash
module load GCC/11.3.0 OpenMPI/4.1.4 ColabFold/1.5.2-CUDA-11.7.0
colabfold_batch --model-type alphafold2_multimer_v3 \
                --num-models 5 --num-recycle 3 <input.fasta> <out>
```

FASTA entries pair the two chains as BINDER_SEQ:TARGET_SEQ on a single line.

Use alphafold2_multimer_v3 throughout. The alphafold2_ptm weights available on
the cluster are corrupted and fail silently rather than raising an error. A100 cards
are recommended for this step.

## Official filter thresholds

Each tool forwards its best designs after its own published filter. For the
AlphaFold2 based chain this mirrors the BindCraft defaults.

```
ipTM           > 0.5
pTM            > 0.55
pLDDT          > 0.8
interface pAE  < 10.85 Å
```

That last figure needs explaining. BindCraft stores the threshold in normalised form
as 0.35, not in ångströms. Undoing the normalisation gives 0.35 multiplied by 31,
which is roughly 10.85 Å. Reading the configuration value as an ångström distance is
an easy mistake to make.

## SLURM notes

```bash
squeue -u $USER
sbatch --array=3,4,17 script.sh
scancel <jobid>
```

Partitions cap at twelve hours, so every long step needs a resume mechanism that
skips outputs already produced. The PyRosetta scoring script writes results line by
line and skips designs already scored, which makes it safe to relaunch after a
timeout. Targets are indexed 1 to 20 and job arrays are the natural way to relaunch
a subset after a partial failure.
