# Cluster 3: Population mismatch and family-structure bias

This folder contains the simulation code for the third bias cluster in the 2SMR bias map manuscript.

## Scientific aim

To show that population mismatch and family-structure bias are structural design problems that are not repaired by generic pleiotropy-robust estimation.

## Main message

Population alignment and within-family designs address different validity threats from those targeted by standard summary-data robustness methods.

## Script

- `run_cluster3_simulation.R`

## Core simulation settings

- Number of replications per scenario: 50
- Scenario parameter grid: 0, 0.25, 0.5, 0.75, 1.0
- Number of SNPs: 500
- True causal effect: 0.3
- Exposure GWAS sample size: 50000
- Outcome GWAS sample size: 50000

## Simulated scenarios

- Population mismatch
- Family-structure bias / dynastic confounding

## Compared methods

- MR-Egger
- MR-RAPS
- IVW_SamePopFilter
- IVW_WithinFamily
- SimSS_IVW
- SimSS_RAPS
- MRMix
- PCMR
- MR_Horse

## Figure mapping

- Panel I: population mismatch
- Panel J: family-structure bias

## Expected outputs

Store generated intermediate results and figure panels in the `outputs/` folder.
Store run logs, warnings, and method-specific diagnostics in the `logs/` folder.
