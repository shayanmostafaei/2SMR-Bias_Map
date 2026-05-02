# Suggested run order

This file summarizes the recommended execution order for the three main simulation domains in the 2SMR Bias Map project.

## Cluster 1
Run:

`simulations/cluster1_overlap_weak_winners/run_cluster1_simulation.R`

Outputs:
- panel-level summaries for overlap, weak instruments, and threshold-based instrument selection
- manuscript-linked figure outputs for Cluster 1

## Cluster 2
Run:

`simulations/cluster2_pleiotropy_ld/run_cluster2_simulation.R`

Outputs:
- panel-level summaries for five pleiotropy / LD architectures
- manuscript-linked figure outputs for Cluster 2
- diagnostic report for method disagreement / concordance

## Cluster 3
Run:

`simulations/cluster3_population_family/run_cluster3_simulation.R`

Outputs:
- panel-level summaries for population mismatch and family-structure bias
- manuscript-linked figure outputs for Cluster 3
- diagnostic report

## Final figure workflow
After all three simulation clusters have been run:

1. check outputs in each cluster folder
2. inspect cluster-specific figure panels
3. assemble the final manuscript-ready merged figure
4. confirm panel labels and captions match the manuscript

## Notes
- Re-run simulations after any substantive method or parameter change
- Keep final manuscript figures separate from temporary development outputs
- Save any manually edited figure assets in a clearly named folder under `figures/`
