# PACE, Protein Affinity Comparative Evaluator

A reproducible HPC pipeline for the method-independent comparison of five de novo
protein binder design tools across 20 therapeutic targets.

MSc Bioinformatics thesis, University of Geneva, 2026. Théo Syrvet.

## The problem

Every de novo binder design tool ships with its own scoring system, and those
systems are not comparable. RFdiffusion scores with AlphaFold2, BindCraft with its
own AlphaFold2 Multimer filters, RFantibody with RoseTTAFold2, BoltzGen with a
composite derived from Boltz. Comparing published pass rates across tools therefore
compares predictors rather than designs.

There is a second problem. Each tool is judged by a model closely related to the one
that generated its designs, and a model tends to be lenient towards its own output.

PACE addresses both by submitting every design, whatever produced it, to a single
independent judge, and then testing that judgement against physics.

## The approach

Evaluation runs on three levels.

**Level 1, official scoring.** Each of the five tools (RFdiffusion, BindCraft,
BoltzGen, PPDiff, RFantibody) generates binders for the same 20 targets, with the
same imposed hotspots, using its own published protocol. Each then filters its own
designs with its own system, and forwards its best 100 per target.

**Level 2, PACE.** Every surviving design is re-evaluated by the same independent
judge, built from two components that ask different questions of the same interface.
OpenFold3 refolds each binder and target complex and reports the chain pair ipTM,
written cpi, which says whether the interface is plausible to a deep learning model.
PyRosetta then scores that same structure on physical grounds, returning the
interface binding free energy dG, the shape complementarity SC and the buried
surface area dSASA, which say whether it is solid. Neither judge is sufficient
alone, and the strength of the method lies in crossing the two.

OpenFold3 is an open reproduction of the AlphaFold3 architecture, trained
separately, and it is not the scoring model of any tool in the benchmark. That
independence is real but partial, since OpenFold3 shares an architectural family
with BoltzGen. Level 3 exists partly to test whether that residual bias matters.

**Level 3, molecular dynamics.** The best design of each target and method pair is
simulated in explicit solvent and judged on physics alone.

### Acceptance criteria

A design must satisfy all three physical thresholds at once. The interface binding
free energy must fall below -15 REU, deliberately stricter than the default of 0 so
that only clearly favourable interfaces survive. Shape complementarity, in the sense
of Lawrence and Colman, must exceed 0.55. Buried interface surface area must exceed
800 Å². Failing any one of the three discards the design, whatever its predicted
confidence.

The thresholds come from the published de novo binder literature. They were not
tuned on this dataset, which keeps them independent of the methods being evaluated.

cpi is deliberately not used as a cutoff. Its distribution is unimodal, with no
natural boundary between a plausible design and one to discard, so a threshold would
cut an arbitrary line through a continuum. It serves for ranking only.

The two measures are weakly correlated by design, which is the point. Each carries
information the other does not, so the combination is used to eliminate poor
designs rather than to separate the best ones.

### Molecular dynamics validation

Static metrics judge a frozen interface. They say nothing about whether it holds
once the system is allowed to move, and a predicted interface can look excellent yet
come apart immediately.

The best design for each target and method therefore runs through all-atom
molecular dynamics in explicit solvent, with the CHARMM36m force field, TIP3P water
and a dodecahedral box neutralised at 150 mM. Each system is simulated in three
independent replicas of 20 ns with randomised initial velocities. Two quantities are
extracted per trajectory, the interface RMSD and the binding free energy estimated
by MM-PBSA, averaged across replicas.

This temporal validation layer is absent from published binder comparison protocols
and is the main methodological contribution of the work.

## What the pipeline showed

Three findings, stated qualitatively while the manuscript is under revision.

**Internal scores are misleading.** The ranking produced by each tool's own metrics
differs substantially from the ranking produced by an independent judge. One method
rates itself very favourably and finishes last under PACE. Another dominates without
anything in its native metrics predicting it.

**The hierarchy is not an artefact of the judge.** Molecular dynamics shares no
architecture with either the generative methods or with PACE, and rests on a force
field rather than on learned weights. It nevertheless recovers the same leading
methods, which rules out the concern that the result was driven by architectural
kinship between OpenFold3 and one of the tools. Interface stability and binding
energy also agree with each other, two independent properties pointing the same way.

**The ranking metric has a scope.** Tested against real biological complexes, cpi
scores natural binders poorly and ranks the true complex first in only three cases
out of eleven. Binder length spans a factor of roughly 85 across that set, and cpi
normalises by residue count, so the metric suits homogeneous de novo binders of
around one hundred residues rather than heterogeneous natural partners. This does
not affect the benchmark, whose designs are of comparable size, but it does bound
where the method can be applied, and it is the reason a third independent level of
validation was added.

## Repository layout

```
env/                    # weight download and environment setup

pipeline/
    01_design/          # input preparation, format conversion, hotspot fixes
    02_sequence/        # FASTA generation, extraction, top N selection
    03_structure/       # post folding structure processing
    04_metrics/         # PyRosetta scoring, SLURM array jobs, one per method
    05_pace/            # aggregation, ranking, common scoring

utils/                  # pipeline health checks, GPU hour accounting
docs/workflow.md        # execution notes and per tool caveats
```

This repository contains the orchestration and evaluation layer only. The five
design tools are third-party software with their own licences and are not vendored
here. Clone them separately, following the workflow notes under docs.

No results or structural data are included, by design.

## Running it

The pipeline targets a SLURM cluster with GPU partitions, running Rocky Linux 9. A
twelve hour wall time limit per job makes resumable scripts a requirement rather
than a convenience, and the 20 targets are indexed 1 to 20 through job arrays.

The tools in the benchmark depend on mutually incompatible software stacks, so each
stage runs in its own conda environment. Exported specifications for all eight sit
under env/conda and are recreated with conda env create. Backbone generation uses
rfdiffusion and SE3nv, sequence design and its AlphaFold2 validation use
af2_binder_design, the remaining generative tools use BindCraft and boltzgen, PACE
refolding uses openfold3, and the dynamics stage uses gmxmmpbsa and mdanalysis.

Two environment variables must be set.

```bash
export PACE_ROOT=/path/to/your/working/directory
export PACE_PYTHON=~/.conda/envs/BindCraft/bin/python
```

Then, from the repository root:

```bash
mkdir -p logs
sbatch pipeline/04_metrics/run_pyrosetta_array.sh
```

SLURM output directives are relative to the submission directory, so the logs
folder must exist before submitting.

Because runs are long and interruptions routine, the utils folder holds a health
check that reports the completion state of each target, and a GPU accounting script
that reconstructs compute time per method and per target from scheduler records.
That accounting is what made the cost side of the comparison a measured quantity
rather than an estimate.

Every tool in the benchmark has undocumented failure modes, and they cost real
debugging time. RFdiffusion needed a CUDA conflict resolved and its IGSO3 cache
handled, PPDiff shipped a corrupted checkpoint, BoltzGen required CIF parsing work
and deadlocked in its data loader. The workflow notes under docs cover these tool by
tool. Mutually incompatible software stacks were isolated in four Apptainer
containers alongside conda environments.

## Stack

Python, R, Bash, SLURM job arrays, Apptainer, Docker, PyRosetta, GROMACS,
OpenFold3, ColabFold, Git.

## Status

Thesis submitted and defended in July 2026, graded 6 out of 6. The manuscript is
under revision for publication, so figures, per method rankings and numerical
results are not published here for the time being.

## Licence

Released under the MIT licence. The third-party design tools retain their own.
