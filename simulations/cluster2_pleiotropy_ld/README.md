# Cluster 2: Pleiotropy and linkage disequilibrium (LD)

This folder contains the simulation code for the second bias cluster in the 2SMR bias map manuscript.

## Scientific aim

To compare estimator behavior across multiple pleiotropy architectures and demonstrate that no single estimator is uniformly best across all pleiotropic settings.

## Main message

Different pleiotropy structures require different responses. Agreement across methods that fail differently is more informative than reliance on a single “robust” estimator.

## Script

- `run_cluster2_simulation.R`

## Core simulation settings

- Number of replications per scenario: 20
- Scenario parameter grid: 0, 0.5, 1.0, 1.5, 2.0
- Number of SNPs: 200
- True causal effect: 0.3
- Exposure GWAS sample size: 50000
- Outcome GWAS sample size: 50000

## Simulated pleiotropy architectures

- Outlier pleiotropy
- Directional pleiotropy
- Balanced weak-instrument pleiotropy
- Correlated pleiotropy / InSIDE violation
- LD-structured pleiotropy

## Figure mapping

- Panel D: outlier pleiotropy
- Panel E: directional pleiotropy
- Panel F: balanced weak-instrument pleiotropy
- Panel G: correlated pleiotropy / InSIDE violation
- Panel H: LD-structured pleiotropy

## Expected outputs

Store generated intermediate results and figure panels in the `outputs/` folder.
Store run logs, warnings, and method-specific diagnostics in the `logs/` folder.
