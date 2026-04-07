# 2SMR Bias Map

This repository contains the simulation code, findings/figures, and planning materials for a methods-oriented workflow on the credibility of two-sample Mendelian randomization (2SMR).

## Project goal

The central claim of this project is that no single Mendelian randomization method can overcome all major bias classes in 2SMR. Instead, credible inference requires **bias-to-method matching**.

We organize the project around major bias clusters:

1. Winner's curse, weak instruments, and participant overlap
2. Pleiotropy and linkage disequilibrium (LD)
3. Population mismatch and family-structure bias
4. Selection bias, index-event bias, and survival bias
5. Interpretation problems, including time-varying and binary exposures

## Repository structure

- `figures/` — final figure panels and figure-specific outputs
- `simulations/` — all simulation scripts organized by bias cluster
- `docs/` — notes, figure plans, and bias-to-method summaries
- `.github/` — issue templates, PR template, and contributor workflow files

## Reproducibility principles

- Each simulation cluster should have:
  - one main script
  - one parameter file
  - one output folder
  - one short README describing the data-generating mechanism
- Every figure should be reproducible from code stored in this repository

## Main workflow

1. Define the bias cluster
2. Specify the causal estimand and failure mode
3. Choose comparison methods
4. Run simulations
5. Create figure panels
6. Record limitations and next-step issues

## Planned simulation roadmap

### Cluster 1
Winner's curse + weak instruments + participant overlap  
Goal: show that overlap correction and weak-instrument robustness solve different parts of the same bias cluster.

### Cluster 2
Pleiotropy + LD  
Goal: show that different pleiotropy structures require different estimators and that no single “robust” method is universally optimal.

### Cluster 3
Population mismatch + family structure  
Goal: show when summary-data 2SMR becomes a design problem rather than a sensitivity-analysis problem.

### Cluster 4
Selection, survival, and index-event bias  
Goal: show how apparently valid instruments can still yield biased causal estimates under non-random sampling.

### Cluster 5
Interpretation and estimand mismatch  
Goal: show when univariable 2SMR no longer matches the scientific question.

## Citation

If you use material from this repository, please cite the associated manuscript once available.

## Contact

Maintainer: Dr. Shayan Mostafaei
Contributors: Ali Taghavi Rad
GitHub: https://github.com/shayanmostafaei     https://github.com/AliTaghavirad 
