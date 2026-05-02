# ============================================================================
# Cluster 1: Winner's curse, weak instruments, and participant overlap
# Methods: Classical IVW, Classical MR-RAPS, MR-SimSS (IVW), MR-SimSS (RAPS)
# ============================================================================
#
# Author: Shayan Mostafaei
# Updated: 2026-04-07
#
# Description:
#   Simulation study for Bias Cluster 1 in the 2SMR Bias Map project.
#   Quantifies bias under participant overlap, weak instruments, and
#   threshold-based instrument selection (winner's curse).
#
# Manuscript mapping:
#   Figure 1A-C
#
# Inputs:
#   - No external inputs required; all data are simulated
#
# Outputs:
#   - cluster1_bias_cluster_figure.pdf
#   - cluster1_bias_cluster_figure.png
#   - cluster1_panelA_summary.csv
#   - cluster1_panelB_summary.csv
#   - cluster1_panelC_summary.csv
#
# Notes:
#   - Uses genome-wide simulations to estimate overlap-driven covariance
#   - Uses MR-SimSS for sample-splitting correction
#   - Uses MR-RAPS for weak-instrument-robust estimation
#   - Reproducible with fixed RNG seed
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
N_REPS <- 50

N_SNPS <- 200000
PROP_EFFECT <- 0.012
N_X <- 200000
N_Y <- 200000

# Repo-aware output paths
base_dir <- if (exists("base_dir")) base_dir else here::here()
out_dir <- if (exists("out_dir")) {
  out_dir
} else {
  here::here("simulations", "cluster1_overlap_weak_winners", "outputs")
}
fig_dir <- if (exists("fig_dir")) {
  fig_dir
} else {
  here::here("figures", "cluster1")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(SEED)

# ============================================================================
# Section 1: Simulation function
# ============================================================================
simulate_gwas_simss <- function(
    n_snps = N_SNPS,
    prop_effect = PROP_EFFECT,
    h2 = 0.5,
    n_x = N_X,
    n_y = N_Y,
    frac_overlap = 0,
    cor_xy = 0.5,
    beta_xy = TRUE_EFFECT
) {
  n_overlap <- frac_overlap * min(n_x, n_y)
  maf <- runif(n_snps, 0.01, 0.5)

  effect_snps <- n_snps * prop_effect
  index <- sample(seq_len(n_snps), ceiling(effect_snps), replace = FALSE)

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
  cov_array[1, 2, ] <- cov_gx_gy
  cov_array[2, 2, ] <- var_gy

  summary_stats <- apply(
    cov_array,
    3,
    function(x) MASS::mvrnorm(n = 1, mu = c(0, 0), Sigma = x)
  )

  beta.exposure <- summary_stats[1, ] + beta_gx
  beta.outcome  <- summary_stats[2, ] + beta_gy

  data.frame(
    SNP = seq_len(n_snps),
    beta.exposure = beta.exposure,
    beta.outcome = beta.outcome,
    se.exposure = sqrt(var_gx),
    se.outcome = sqrt(var_gy)
  )
}

# ============================================================================
# Section 2: Estimator application
# ============================================================================
apply_estimators <- function(data, f_threshold = 10, overlap_frac = 0) {
  pval_thresh <- pchisq(f_threshold, df = 1, lower.tail = FALSE)
  num_splits <- ifelse(overlap_frac == 0, 2, 3)

  data_sig <- data |>
    filter(2 * pnorm(abs(beta.exposure / se.exposure), lower.tail = FALSE) < pval_thresh)

  # Classical IVW
  if (nrow(data_sig) > 2) {
    res_ivw <- summary(
      lm(beta.outcome ~ -1 + beta.exposure, weights = 1 / se.outcome^2, data = data_sig)
    )
    b_class_ivw <- res_ivw$coef[1, 1]
  } else {
    b_class_ivw <- NA_real_
  }

  # Classical MR-RAPS
  if (nrow(data_sig) > 2) {
    suppressWarnings(suppressMessages({
      res_raps <- mr.raps(
        data_sig$beta.exposure,
        data_sig$beta.outcome,
        data_sig$se.exposure,
        data_sig$se.outcome,
        over.dispersion = FALSE
      )
      b_class_raps <- res_raps$beta.hat
    }))
  } else {
    b_class_raps <- NA_real_
  }

  # MR-SimSS (IVW)
  simss_ivw <- mr_simss(
    data,
    mr_method = "mr_ivw",
    threshold = pval_thresh,
    splits = num_splits,
    n.iter = 1000
  )
  b_simss_ivw <- simss_ivw$summary$b

  # MR-SimSS (RAPS)
  simss_raps <- mr_simss(
    data,
    mr_method = "mr_raps",
    threshold = pval_thresh,
    splits = num_splits,
    n.iter = 1000
  )
  b_simss_raps <- simss_raps$summary$b

  c(
    Classical_IVW  = b_class_ivw,
    Classical_RAPS = b_class_raps,
    SimSS_IVW      = b_simss_ivw,
    SimSS_RAPS     = b_simss_raps
  )
}

# ============================================================================
# Section 3: Run simulations
# ============================================================================
message("Running Cluster 1 simulations ...")

# Panel A: participant overlap
message("Panel A: Participant overlap")
overlap_grid <- seq(0, 1, by = 0.25)
res_A_list <- vector("list", length(overlap_grid) * N_REPS)
idx <- 1

for (ov in overlap_grid) {
  message("  overlap = ", ov)
  for (i in seq_len(N_REPS)) {
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

# Panel B: weak instruments
message("Panel B: Weak instruments under complete overlap")
h2_grid <- seq(0.05, 0.5, length.out = 5)
res_B_list <- vector("list", length(h2_grid) * N_REPS)
idx <- 1

for (h2_val in h2_grid) {
  message("  h2 = ", round(h2_val, 3))
  for (i in seq_len(N_REPS)) {
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

# Panel C: winner's curse / threshold-based selection
message("Panel C: Threshold-based instrument selection")
threshold_grid <- c(5, 10, 20, 30, 40)
res_C_list <- vector("list", length(threshold_grid) * N_REPS)
idx <- 1

for (thresh in threshold_grid) {
  message("  F-threshold = ", thresh)
  for (i in seq_len(N_REPS)) {
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
# Section 4: Summaries for plotting and export
# ============================================================================
make_summary <- function(data) {
  data |>
    group_by(Param, Method) |>
    summarise(
      Mean = mean(Estimate, na.rm = TRUE),
      SD = sd(Estimate, na.rm = TRUE),
      SE = SD / sqrt(sum(!is.na(Estimate))),
      N = sum(!is.na(Estimate)),
      .groups = "drop"
    ) |>
    na.omit()
}

sum_A <- make_summary(res_A)
sum_B <- make_summary(res_B)
sum_C <- make_summary(res_C)

write.csv(sum_A, file.path(out_dir, "cluster1_panelA_summary.csv"), row.names = FALSE)
write.csv(sum_B, file.path(out_dir, "cluster1_panelB_summary.csv"), row.names = FALSE)
write.csv(sum_C, file.path(out_dir, "cluster1_panelC_summary.csv"), row.names = FALSE)

# ============================================================================
# Section 5: Plotting
# ============================================================================
message("Generating figure ...")

# Colorblind-safe palette
my_colors <- c(
  Classical_IVW  = "#000000",
  Classical_RAPS = "#0072B2",
  SimSS_IVW      = "#D55E00",
  SimSS_RAPS     = "#009E73"
)

method_labels <- c(
  Classical_IVW  = "Classical IVW",
  Classical_RAPS = "Classical MR-RAPS",
  SimSS_IVW      = "MR-SimSS (IVW)",
  SimSS_RAPS     = "MR-SimSS (RAPS)"
)

method_linetypes <- c(
  Classical_IVW  = "solid",
  Classical_RAPS = "dashed",
  SimSS_IVW      = "dotted",
  SimSS_RAPS     = "solid"
)

method_shapes <- c(
  Classical_IVW  = 16,
  Classical_RAPS = 17,
  SimSS_IVW      = 15,
  SimSS_RAPS     = 18
)

plot_panel <- function(df_summary, x_label, title, invert_x = FALSE) {
  p <- ggplot(
    df_summary,
    aes(
      x = Param,
      y = Mean,
      color = Method,
      linetype = Method,
      shape = Method,
      group = Method
    )
  ) +
    geom_hline(
      yintercept = TRUE_EFFECT,
      linetype = "dashed",
      color = "grey25",
      linewidth = 0.8
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(ymin = Mean - 1.96 * SE, ymax = Mean + 1.96 * SE),
      width = 0.04 * diff(range(df_summary$Param)),
      linewidth = 0.5
    ) +
    scale_color_manual(values = my_colors, labels = method_labels) +
    scale_linetype_manual(values = method_linetypes, labels = method_labels) +
    scale_shape_manual(values = method_shapes, labels = method_labels) +
    labs(
      title = title,
      x = x_label,
      y = "Estimated causal effect"
    ) +
    theme_classic(base_family = "sans", base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "none"
    )

  if (invert_x) {
    p <- p + scale_x_reverse()
  }

  p
}

pA <- plot_panel(
  sum_A,
  x_label = "Sample overlap (%)",
  title = "A. Participant overlap"
)

pB <- plot_panel(
  sum_B,
  x_label = "Exposure heritability (proxy for instrument strength)",
  title = "B. Weak instruments under complete overlap",
  invert_x = TRUE
)

pC <- plot_panel(
  sum_C,
  x_label = "Selection threshold (F-statistic)",
  title = "C. Threshold-based instrument selection"
)

legend_plot <- ggplot(
  sum_A,
  aes(
    x = Param,
    y = Mean,
    color = Method,
    linetype = Method,
    shape = Method
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = my_colors, labels = method_labels) +
  scale_linetype_manual(values = method_linetypes, labels = method_labels) +
  scale_shape_manual(values = method_shapes, labels = method_labels) +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  )

shared_legend <- get_legend(legend_plot)

top_row <- plot_grid(pA, pB, pC, ncol = 3, align = "hv")
final_plot <- plot_grid(top_row, shared_legend, ncol = 1, rel_heights = c(1, 0.12))

# ============================================================================
# Section 6: Save outputs
# ============================================================================
pdf_name <- "cluster1_bias_cluster_figure.pdf"
png_name <- "cluster1_bias_cluster_figure.png"

ggsave(
  filename = file.path(fig_dir, pdf_name),
  plot = final_plot,
  width = 16,
  height = 6,
  device = cairo_pdf
)

ggsave(
  filename = file.path(fig_dir, png_name),
  plot = final_plot,
  width = 16,
  height = 6,
  dpi = 500
)

message(
  "\n", strrep("=", 72), "\n",
  "Cluster 1 complete [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]\n",
  strrep("-", 72), "\n",
  "Summary outputs : ", out_dir, "\n",
  "Figure outputs  : ", fig_dir, "\n",
  "Files written   :\n",
  "  ", file.path(out_dir, "cluster1_panelA_summary.csv"), "\n",
  "  ", file.path(out_dir, "cluster1_panelB_summary.csv"), "\n",
  "  ", file.path(out_dir, "cluster1_panelC_summary.csv"), "\n",
  "  ", file.path(fig_dir, pdf_name), "\n",
  "  ", file.path(fig_dir, png_name), "\n",
  strrep("=", 72), "\n"
)
