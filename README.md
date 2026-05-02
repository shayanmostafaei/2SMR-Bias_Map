# 2SMR Bias Map

This repository contains the simulation code, figures, and supporting materials for the project **“From scalable to credible: a practical bias map for two-sample Mendelian randomization.”** The project develops a bias-aware framework for thinking about the credibility of two-sample Mendelian randomization (2SMR), with an emphasis on simulation-based illustration of major validity threats.

## Project goal

The central claim of this project is that **no single Mendelian randomization method is uniformly best across the major bias domains in 2SMR**. Credible inference instead requires **bias-to-method matching**: different methods solve different problems, rely on different assumptions, and answer different causal questions.

The project is organized around three major simulated bias domains used in the current manuscript:

1. **Winner’s curse, weak instruments, and participant overlap**
2. **Pleiotropy and linkage disequilibrium (LD)**
3. **Population mismatch and family-structure bias**

Additional repository materials cover broader design and interpretation issues relevant to 2SMR, including selection bias, index-event bias, and estimand mismatch.

## Manuscript status

This repository supports the manuscript currently titled:

**From scalable to credible: a practical bias map for two-sample Mendelian randomization**

The manuscript is structured as a commentary-style synthesis supported by simulation-based illustrations. The key takeaway is:

> No single MR method is best. Different methods answer different questions, and the appropriate choice depends on the dominant bias structure, the data-generating setting, and the causal estimand.

## Figure map

The current main simulation figure is organized as follows:

- **Panels A–C:** participant overlap, weak instruments, and threshold-based instrument selection
- **Panels D–H:** pleiotropy and LD across five architectures  
  - outlier pleiotropy  
  - directional pleiotropy  
  - balanced weak-instrument pleiotropy  
  - correlated pleiotropy / InSIDE violation  
  - LD-structured pleiotropy
- **Panels I–J:** population mismatch and family-structure confounding

Together, these panels illustrate that overlap correction, weak-instrument robustness, pleiotropy-aware estimation, population matching, and within-family designs address different threats to validity rather than converging on a single universally robust solution.

## Repository structure

- `figures/` — final figure panels, assembled figures, and figure-specific outputs
- `simulations/` — simulation scripts organized by bias domain
- `.github/` — issue templates, pull request templates, and contributor workflow files

## Reproducibility principles

Each simulation domain should be reproducible from code stored in this repository. Where possible, each cluster should include:

- one main script
- one parameter definition or settings file
- one output folder
- one short README describing the data-generating mechanism and estimand

Figures included in the manuscript should be traceable to scripts and saved outputs in this repository.

## Main workflow

1. Define the bias domain
2. Specify the estimand and dominant failure mode
3. Choose representative comparison methods
4. Run the simulations
5. Generate figure panels
6. Record assumptions, limitations, and interpretation notes


## Contact

**Maintainer:** Shayan Mostafaei  
**Contributor:** Ali Taghavi Rad
