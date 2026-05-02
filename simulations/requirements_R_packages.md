# Required R packages

The simulation scripts in this repository use some or all of the following R packages.

## Core packages
- MASS
- dplyr
- tidyr
- ggplot2
- cowplot
- Matrix
- here

## Mendelian randomization / method-specific packages
- MendelianRandomization
- mr.raps
- mr.simss
- mr.divw
- MRMix
- RBMR
- MR.LDP
- cause
- PCMR
- MRSamePopTest

## Parallel / workflow support
- doParallel
- foreach

## Bayesian / external model support
- R2jags

## External software notes
Some scripts may require:
- **JAGS** installed and accessible from R
- successful installation of method-specific packages that may not be available on standard CRAN workflows

## Recommended setup practice
For reproducibility, consider using:
- `renv`
- a fixed R version
- a saved session info log after each major simulation run

## Suggested verification
Before running the full simulation pipeline, test package loading interactively in R with:

```r
packages <- c(
  "MASS", "dplyr", "tidyr", "ggplot2", "cowplot", "Matrix", "here",
  "MendelianRandomization", "mr.raps", "mr.simss", "mr.divw",
  "MRMix", "RBMR", "MR.LDP", "cause", "PCMR", "MRSamePopTest",
  "doParallel", "foreach", "R2jags"
)

sapply(packages, require, character.only = TRUE)
