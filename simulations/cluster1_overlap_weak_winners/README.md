# Cluster 1: Winner’s curse, weak instruments, and participant overlap

This folder contains the simulation code for the first bias cluster in the 2SMR bias map manuscript.

## Scientific aim

To illustrate that winner’s curse, weak-instrument bias, and participant overlap form an interacting bias cluster rather than three separable caveats.

## Main message

No single correction solves this cluster on its own. Simulated sample splitting addresses overlap-related bias, but weak-instrument attenuation remains unless paired with a weak-instrument-robust estimator.

## Script

- `run_cluster1_simulation.R`

## Core simulation settings

- Number of repetitions per panel: 50
- Number of SNPs: 200000
- True causal effect: 0.3
- Exposure GWAS sample size: 200000
- Outcome GWAS sample size: 200000
- Default heritability: 0.5

## Compared methods

- Classical IVW
- Classical MR-RAPS
- MR-SimSS (IVW)
- MR-SimSS (RAPS)

## Figure mapping

- Panel A: participant overlap
- Panel B: weak instruments
- Panel C: threshold-based instrument selection / winner’s curse

## Expected outputs

Store generated intermediate results and figure panels in the `outputs/` folder.
Store run logs or warnings in the `logs/` folder.
