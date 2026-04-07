# ============================================================================
# Deconstructing the "Bias Cluster 1"
# Methods: Classical IVW, Classical MR-RAPS, SimSS-IVW, SimSS-RAPS
# ============================================================================
#
# Author: Shayan Mostafaei
# Date: 2026-04-07
#
# Description:
#   Simulation study to quantify bias under sample overlap, weak instruments,
#   and winner's curse. Produces a multi-panel figure comparing classical and
#   SimSS estimators.
#
# Inputs:
#   - No external inputs required. All data are simulated.
#
# Outputs:
#   - Fig1_Bias_Cluster_HighRes.pdf
#   - Fig1_Bias_Cluster_HighRes.png
#
# Notes:
#   - Uses genome-wide simulations to estimate lambda and instrument strength.
#   - Overlap, weak IVs, and selection thresholds are varied across panels.
#   - Results are reproducible given the fixed RNG seed.
# ============================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(mr.raps)
  library(mr.simss)
  library(here)
})

# ============================================================================
# Section 0: Global settings
# ============================================================================
SEED <- 20260406
TRUE_EFFECT <- 0.3
N_REPS <- 100

# Allow users to override out_dir before sourcing this script
out_dir <- if (exists("out_dir")) out_dir else here::here("outputs", "simulations")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

# ============================================================================
# Section 1: Simulation function
# ============================================================================
#' Simulate GWAS summary statistics for MR-SimSS
#'
#' @param n_snps Number of SNPs to simulate.
#' @param prop_effect Proportion of SNPs with non-zero effects.
#' @param h2 Heritability (instrument strength).
#' @param n_x Sample size for exposure GWAS.
#' @param n_y Sample size for outcome GWAS.
#' @param frac_overlap Fraction of sample overlap.
#' @param cor_xy Observational confounding/correlation.
#' @param beta_xy True causal effect.
#' @return A data frame with SNP, beta, and standard error columns.
#'
#' @noRd
simulate_gwas_simss <- function(
    n_snps = 1000000,
    prop_effect = 0.012,
    h2 = 0.5,
    n_x = 200000,
    n_y = 200000,
    frac_overlap = 0,
    cor_xy = 0.5,
    beta_xy = 0.3
) {
  n_overlap <- frac_overlap * min(n_x, n_y)
  maf <- runif(n_snps, 0.01, 0.5)

  effect_snps <- n_snps * prop_effect
  index <- sample(1:n_snps, ceiling(effect_snps), replace = FALSE)

  beta_gx <- rep(0, n_snps)
  beta_gx[index] <- rnorm(length(index), 0, 1)

  var_x <- sum(2 * maf * (1 - maf) * beta_gx^2) / h2
  if (var_x != 0) {
    beta_gx <- beta_gx / sqrt(var_x)
  }

  beta_gy <- beta_gx * beta_xy

  var_gx <- 1 / (n_x * 2 * maf * (1 - maf))
  var_gy <- 1 / (n_y * 2 * maf * (1 - maf))
  cov_gx_gy <- ((n_overlap * cor_xy) / (n_x * n_y)) * (1 / (2 * maf * (1 - maf)))

  cov_array <- array(dim = c(2, 2, n_snps))
  cov_array[1, 1, ] <- var_gx
  cov_array[2, 1, ] <- cov_gx_gy
  cov_array[1, 2, ] <- cov_array[2, 1, ]
  cov_array[2, 2, ] <- var_gy

  summary_stats <- apply(
    cov_array,
    3,
    function(x) MASS::mvrnorm(n = 1, mu = c(0, 0), Sigma = x)
  )

  beta.exposure <- summary_stats[1, ] + beta_gx
  beta.outcome <- summary_stats[2, ] + beta_gy

  data.frame(
    SNP = 1:n_snps,
    beta.exposure = beta.exposure,
    beta.outcome = beta.outcome,
    se.exposure = sqrt(var_gx),
    se.outcome = sqrt(var_gy)
  )
}

# ============================================================================
# Section 2: Estimator application
# ============================================================================
#' Apply classical and SimSS estimators
#'
#' @param data GWAS summary statistics from simulate_gwas_simss().
#' @param f_threshold F-statistic threshold used for SNP selection.
#' @param overlap_frac Sample overlap fraction (controls split strategy).
#' @return Named numeric vector of estimates.
#'
#' @noRd
apply_estimators <- function(data, f_threshold = 10, overlap_frac = 0) {
  pval_thresh <- pchisq(f_threshold, df = 1, lower.tail = FALSE)
  num_splits <- ifelse(overlap_frac == 0, 2, 3)

  data_sig <- data |> 
    filter(2 * (pnorm(abs(beta.exposure / se.exposure), lower.tail = FALSE)) < pval_thresh)

  if (nrow(data_sig) > 2) {
    res.ivw <- summary(
      lm(beta.outcome ~ -1 + beta.exposure, weights = 1 / se.outcome^2, data = data_sig)
    )
    b_class_ivw <- res.ivw$coef[1, 1]
  } else {
    b_class_ivw <- NA_real_
  }

  if (nrow(data_sig) > 2) {
    suppressWarnings(suppressMessages({
      res.raps <- mr.raps(
        data_sig$beta.exposure,
        data_sig$beta.outcome,
        data_sig$se.exposure,
        data_sig$se.outcome,
        over.dispersion = FALSE
      )
      b_class_raps <- res.raps$beta.hat
    }))
  } else {
    b_class_raps <- NA_real_
  }

  simss_ivw <- mr_simss(
    data,
    mr_method = "mr_ivw",
    threshold = pval_thresh,
    splits = num_splits,
    n.iter = 1000
  )
  b_simss_ivw <- simss_ivw$summary$b

  simss_raps <- mr_simss(
    data,
    mr_method = "mr_raps",
    threshold = pval_thresh,
    splits = num_splits,
    n.iter = 1000
  )
  b_simss_raps <- simss_raps$summary$b

  c(
    Classical_IVW = b_class_ivw,
    Classical_RAPS = b_class_raps,
    SimSS_IVW = b_simss_ivw,
    SimSS_RAPS = b_simss_raps
  )
}

# ============================================================================
# Section 3: Simulation execution
# ============================================================================
message("Running Panel A: Participant Overlap ...")
res_A_list <- vector("list", length(seq(0, 1, by = 0.25)) * N_REPS)
idx <- 1
for (ov in seq(0, 1, by = 0.25)) {
  message("  Overlap: ", ov)
  for (i in 1:N_REPS) {
    d <- simulate_gwas_simss(frac_overlap = ov, h2 = 0.5, beta_xy = TRUE_EFFECT)
    est <- apply_estimators(d, f_threshold = 10, overlap_frac = ov)
    res_A_list[[idx]] <- data.frame(
      Param = ov * 100,
      Method = names(est),
      Estimate = est
    )
    idx <- idx + 1
  }
}
res_A <- bind_rows(res_A_list)

message("Running Panel B: Weak Instruments ...")
res_B_list <- vector("list", length(seq(0.05, 0.5, length.out = 5)) * N_REPS)
idx <- 1
for (h2_val in seq(0.05, 0.5, length.out = 5)) {
  message("  Heritability: ", h2_val)
  for (i in 1:N_REPS) {
    d <- simulate_gwas_simss(frac_overlap = 1.0, h2 = h2_val, beta_xy = TRUE_EFFECT)
    est <- apply_estimators(d, f_threshold = 10, overlap_frac = 1.0)
    res_B_list[[idx]] <- data.frame(
      Param = h2_val,
      Method = names(est),
      Estimate = est
    )
    idx <- idx + 1
  }
}
res_B <- bind_rows(res_B_list)

message("Running Panel C: Winner's Curse ...")
res_C_list <- vector("list", length(c(5, 10, 20, 30, 40)) * N_REPS)
idx <- 1
for (thresh in c(5, 10, 20, 30, 40)) {
  message("  F-Threshold: ", thresh)
  for (i in 1:N_REPS) {
    d <- simulate_gwas_simss(frac_overlap = 0, h2 = 0.3, beta_xy = TRUE_EFFECT)
    est <- apply_estimators(d, f_threshold = thresh, overlap_frac = 0)
    res_C_list[[idx]] <- data.frame(
      Param = thresh,
      Method = names(est),
      Estimate = est
    )
    idx <- idx + 1
  }
}
res_C <- bind_rows(res_C_list)

# ============================================================================
# Section 4: Plotting
# ============================================================================
message("Generating high-resolution plots ...")

my_colors <- c(
  Classical_IVW = "#D55E00",
  Classical_RAPS = "#CC79A7",
  SimSS_IVW = "#E69F00",
  SimSS_RAPS = "#0072B2"
)

method_labels <- c(
  Classical_IVW = "1. Classical IVW",
  Classical_RAPS = "2. Classical MR-RAPS",
  SimSS_IVW = "3. MR-SimSS (IVW)",
  SimSS_RAPS = "4. MR-SimSS (RAPS)"
)

#' Build a panel plot for a simulation scenario
#'
#' @param data Simulation results data frame.
#' @param x_label Label for the x-axis.
#' @param title Plot title.
#' @param is_x_discrete Whether the x-axis should be treated as discrete.
#' @param invert_x Whether to reverse the x-axis.
#' @return A ggplot object.
#'
#' @noRd
plot_panel <- function(data, x_label, title, is_x_discrete = FALSE, invert_x = FALSE) {
  df_summary <- data |
    group_by(Param, Method) |
    summarise(
      Mean = mean(Estimate, na.rm = TRUE),
      SE = sd(Estimate, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) |
    na.omit()

  err_width <- if (is_x_discrete) 0.5 else 0.05 * max(df_summary$Param)

  p <- ggplot(df_summary, aes(x = Param, y = Mean, color = Method, group = Method)) +
    geom_hline(
      yintercept = TRUE_EFFECT,
      linetype = "dashed",
      color = "darkgreen",
      linewidth = 1.2
    ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = Mean - 1.96 * SE, ymax = Mean + 1.96 * SE),
      width = err_width,
      linewidth = 0.8
    ) +
    scale_color_manual(values = my_colors, labels = method_labels) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12)
    ) +
    labs(title = title, y = "Estimated Causal Effect", x = x_label)

  if (invert_x) {
    p <- p + scale_x_reverse()
  }

  p
}

pA <- plot_panel(res_A, "Sample Overlap (%)", "A) Participant Overlap Bias")

pB <- plot_panel(
  res_B,
  "Instrument Strength (Heritability)",
  "B) Weak Instrument Bias\n(Condition: 100% Overlap)",
  invert_x = TRUE
)

pC <- plot_panel(
  res_C,
  "Selection Threshold (F-statistic)",
  "C) Winner's Curse\n(Condition: 0% Overlap)"
)

legend_plot <- ggplot(res_A, aes(x = Param, y = Estimate, color = Method)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = my_colors, labels = method_labels) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 13)
  )

shared_legend <- get_legend(legend_plot)

top_row <- plot_grid(pA, pB, pC, ncol = 3, labels = NULL)
final_plot <- plot_grid(top_row, shared_legend, ncol = 1, rel_heights = c(1, 0.1))

# ============================================================================
# Section 5: Save outputs
# ============================================================================
pdf_name <- "Fig1_Bias_Cluster_HighRes.pdf"
png_name <- "Fig1_Bias_Cluster_HighRes.png"

suppressWarnings({
  ggsave(file.path(out_dir, pdf_name), final_plot, width = 17, height = 6.5, device = cairo_pdf)
  ggsave(file.path(out_dir, png_name), final_plot, width = 17, height = 6.5, dpi = 300)
})

message(
  "\n", strrep("=", 70), "\n",
  "Bias Cluster 1 complete  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]\n",
  strrep("-", 70), "\n",
  "Output directory : ", out_dir, "\n",
  "Files written    :\n",
  "  ", pdf_name, "\n",
  "  ", png_name, "\n",
  strrep("=", 70)
)