# ============================================================================
# Cluster 2: Pleiotropy and linkage disequilibrium (LD)
# Methods:
#   MR-Egger, dIVW, MRMix, MR-SimSS (IVW), MR-SimSS (RAPS), MR-RAPS,
#   RBMR, MR-LDP, PCMR, CAUSE, MR-Horse
# ============================================================================
#
# Author: Shayan Mostafaei
# Updated: 2026-04-07
#
# Description:
#   Simulation study for Bias Cluster 2 in the 2SMR Bias Map project.
#   Compares estimator behavior across five pleiotropy architectures:
#   outlier, directional, balanced weak-instrument, correlated/InSIDE-
#   violating, and LD-structured pleiotropy.
#
# Manuscript mapping:
#   Figure 1D-H
#
# Inputs:
#   - No external data required; all data are simulated
#
# Outputs:
#   - cluster2_pleiotropy_ld_figure.pdf
#   - cluster2_pleiotropy_ld_figure.png
#   - cluster2_panelD_summary.csv
#   - cluster2_panelE_summary.csv
#   - cluster2_panelF_summary.csv
#   - cluster2_panelG_summary.csv
#   - cluster2_panelH_summary.csv
#   - cluster2_diagnostic_report.txt
#
# Notes:
#   - Computationally intensive; several methods may fail in some settings
#   - Includes parallel execution
#   - Saves panel-level summaries and concordance diagnostics
# ============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(MASS)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(MendelianRandomization)
  library(MRMix)
  library(RBMR)
  library(MR.LDP)
  library(mr.divw)
  library(mr.simss)
  library(Matrix)
  library(cause)
  library(doParallel)
  library(foreach)
  library(R2jags)
  library(PCMR)
  library(mr.raps)
  library(here)
})

# ============================================================================
# Section 0: Global configuration
# ============================================================================
SEED_NUMBER <- 20260409

# Simulation scale
N_SIMULATIONS   <- 20
SCENARIO_PARAMS <- seq(0, 2, by = 0.5)

# Data-generation parameters
N_SNPS         <- 200
TRUE_THETA     <- 0.3
SAMPLE_SIZE_N1 <- 50000
SAMPLE_SIZE_N2 <- 50000
BASE_VARIANCE  <- 0.05

# Algorithm / hardware settings
NUM_CORES    <- 4
MAX_ITER_OPT <- 10000

# MR-Horse settings
JAGS_ITER   <- 12000
JAGS_BURNIN <- 5000

# Repo-aware outputs
out_dir <- if (exists("out_dir")) {
  out_dir
} else {
  here::here("simulations", "cluster2_pleiotropy_ld", "outputs")
}
fig_dir <- if (exists("fig_dir")) {
  fig_dir
} else {
  here::here("figures", "cluster2")
}
log_dir <- if (exists("log_dir")) {
  log_dir
} else {
  here::here("simulations", "cluster2_pleiotropy_ld", "logs")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(SEED_NUMBER)

# ============================================================================
# Section 1: Parallel backend
# ============================================================================
if (exists("cl")) {
  try(stopCluster(cl), silent = TRUE)
}
cl <- makeCluster(NUM_CORES)
registerDoParallel(cl)

# ============================================================================
# Section 2: Data generation model
# ============================================================================
generate_pleiotropy_data <- function(param, scenario = "A", n_snps = N_SNPS) {

  block_size <- 20
  num_blocks <- ceiling(n_snps / block_size)
  rho_base <- min(0.65, 0.12 + param * 0.05)

  single_block <- toeplitz(rho_base^(0:(block_size - 1)))
  ld_list <- replicate(num_blocks, single_block, simplify = FALSE)
  ld_matrix <- as.matrix(Matrix::bdiag(ld_list))
  ld_matrix <- ld_matrix[1:n_snps, 1:n_snps]

  sd_alpha <- sqrt(BASE_VARIANCE)
  alpha <- rnorm(n_snps, 0, sd_alpha)
  gamma <- rep(0, n_snps)

  if (scenario == "A") {
    # Outlier pleiotropy
    gamma <- rnorm(n_snps, 0, 0.01)
    outlier_idx <- sample(seq_len(n_snps), size = floor(n_snps * 0.1))
    gamma[outlier_idx] <- rt(length(outlier_idx), df = 3) * (param * 0.3)

  } else if (scenario == "B") {
    # Directional pleiotropy
    bias_drift <- rnorm(1, mean = param * 0.04, sd = 0.01)
    gamma <- rnorm(n_snps, mean = bias_drift, sd = 0.01)

  } else if (scenario == "C") {
    # Balanced weak-instrument pleiotropy
    alpha <- alpha * 0.6
    gamma <- rnorm(n_snps, 0, param * 0.02)

  } else if (scenario == "D") {
    # Correlated pleiotropy / InSIDE violation
    invalid_idx <- sample(seq_len(n_snps), size = floor(n_snps * 0.3))
    gamma[invalid_idx] <- alpha[invalid_idx] * (param * 0.5)
    gamma[-invalid_idx] <- rnorm(n_snps - length(invalid_idx), 0, 0.01)

  } else if (scenario == "E") {
    # LD-structured pleiotropy
    noise_raw <- rnorm(n_snps, 0, 0.06)
    gamma <- as.numeric(ld_matrix %*% noise_raw) * (param * 0.9)
  }

  chol_L <- t(chol(ld_matrix))
  alpha_ld <- as.numeric(chol_L %*% alpha)
  gamma_ld <- as.numeric(chol_L %*% gamma)

  overlap_ratio <- min(0.5, 0.1 + param * 0.1)
  se_corr <- overlap_ratio / sqrt(SAMPLE_SIZE_N1 * SAMPLE_SIZE_N2)
  cov_matrix <- matrix(c(1 / SAMPLE_SIZE_N1, se_corr,
                         se_corr, 1 / SAMPLE_SIZE_N2), 2, 2)
  cov_matrix <- cov_matrix + diag(1e-9, 2)

  beta_exp <- numeric(n_snps)
  beta_out <- numeric(n_snps)

  for (j in seq_len(n_snps)) {
    mu <- c(alpha_ld[j], alpha_ld[j] * TRUE_THETA + gamma_ld[j])
    draw <- MASS::mvrnorm(1, mu = mu, Sigma = cov_matrix)
    beta_exp[j] <- draw[1]
    beta_out[j] <- draw[2]
  }

  maf <- rbeta(n_snps, 1, 3) * 0.4 + 0.05
  se_exp <- sqrt(1 / (2 * maf * (1 - maf) * SAMPLE_SIZE_N1))
  se_out <- sqrt(1 / (2 * maf * (1 - maf) * SAMPLE_SIZE_N2))
  f_stats <- (beta_exp^2) / (se_exp^2)

  dat <- data.frame(
    SNP = paste0("rs", seq_len(n_snps)),
    beta.exposure = beta_exp,
    se.exposure = se_exp,
    beta.outcome = beta_out,
    se.outcome = se_out
  )

  list(dat = dat, ld_matrix = ld_matrix, mean_f = mean(f_stats))
}

# ============================================================================
# Section 3: Helper models and patches
# ============================================================================
mr_horse_model_code <- "model {
  for (i in 1:N){
    by[i] ~ dnorm(mu[i], 1/(sy[i] * sy[i]))
    mu[i] <- theta * bx0[i] + alpha[i]
    bx[i] ~ dnorm(bx0[i], 1 / (sx[i] * sx[i]))
    bx0[i] ~ dnorm(mx0 + (sqrt(vx0)/(tau * phi[i])) * rho[i] * alpha[i], 1 / ((1 - rho[i]^2) * vx0))
    rho[i] ~ dunif(-1, 1)
    alpha[i] ~ dnorm(0, 1/(tau * tau * phi[i] * phi[i]))
    phi[i] <- a[i] / sqrt(b[i])
    a[i] ~ dnorm(0,1) T(0,)
    b[i] ~ dgamma(0.5,0.5)
  }
  c ~ dnorm(0,1) T(0,)
  d ~ dgamma(0.5,0.5)
  tau <- c / sqrt(d)
  vx0 ~ dnorm(0,1) T(0,)
  mx0 ~ dnorm(0,1)
  theta ~ dunif(-10, 10)
}"

patch_RBMR_func <- function() {
  fn <- RBMR::RBMR_func
  fn_body <- deparse(body(fn))
  fn_body <- c("{", "  mu <- gamma", "  muA <- alpha", fn_body[-1])
  body(fn) <- as.call(parse(text = fn_body))
  fn
}
RBMR_func_patched <- patch_RBMR_func()

METHOD_NAMES <- c(
  "MR_Egger", "dIVW", "MRMix", "SimSS_IVW", "SimSS_RAPS", "MR_RAPS",
  "RBMR", "MR_LDP", "PCMR", "CAUSE", "MR_Horse"
)

# ============================================================================
# Section 4: Core estimation wrapper
# ============================================================================
apply_all_methods <- function(dat, ld_matrix, sim_id) {
  res <- list()
  dat <- dat[complete.cases(dat), ]
  if (nrow(dat) < 10) return(res)

  mr_obj <- tryCatch({
    mr_input(
      bx = dat$beta.exposure,
      bxse = dat$se.exposure,
      by = dat$beta.outcome,
      byse = dat$se.outcome
    )
  }, error = function(e) NULL)
  if (is.null(mr_obj)) return(res)

  res$MR_Egger <- tryCatch(mr_egger(mr_obj)$Estimate, error = function(e) NA)

  res$dIVW <- tryCatch({
    w <- 1 / (dat$se.outcome^2)
    num <- sum(w * dat$beta.exposure * dat$beta.outcome)
    den <- sum(w * (dat$beta.exposure^2 - dat$se.exposure^2))
    num / den
  }, error = function(e) NA)

  res$MRMix <- tryCatch(
    MRMix(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)$theta,
    error = function(e) NA
  )

  res$SimSS_IVW <- tryCatch(
    mr_simss(data = dat, mr_method = "mr_ivw", splits = 3)$summary$b,
    error = function(e) NA
  )

  res$SimSS_RAPS <- tryCatch(
    mr_simss(data = dat, mr_method = "mr_raps", splits = 3)$summary$b,
    error = function(e) NA
  )

  res$MR_RAPS <- tryCatch(
    mr.raps(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)$beta.hat,
    error = function(e) NA
  )

  res$RBMR <- tryCatch({
    p <- nrow(dat)
    RBMR_func_patched(
      bh1 = dat$beta.exposure,
      bh2 = dat$beta.outcome,
      se1 = dat$se.exposure,
      se2 = dat$se.outcome,
      R = ld_matrix,
      gamma = rep(0.01, p),
      alpha = rep(0.01, p),
      sgga2 = 0.01,
      sgal2 = 0.01,
      beta0 = 0,
      alphag = 8,
      betag = 4,
      constr = 0,
      epsStopLogLik = 1e-7,
      maxIter = MAX_ITER_OPT
    )$beta0
  }, error = function(e) NA)

  res$MR_LDP <- tryCatch({
    p <- nrow(dat)
    invisible(capture.output(
      mrldp_res <- MRLDP_SimPXvb(
        dat$beta.exposure, dat$beta.outcome,
        dat$se.exposure, dat$se.outcome,
        rep(0, p), rep(0, p), 0, 0.01, 0.01,
        ld_matrix,
        constr = 0,
        epsStopLogLik = 1e-7,
        maxIter = MAX_ITER_OPT,
        model = 2
      )
    ))
    mrldp_res$beta0
  }, error = function(e) NA)

  res$PCMR <- tryCatch({
    suppressWarnings(capture.output(
      pcmr_res <- PCMR(
        beta_ex = dat$beta.exposure,
        beta_ex_se = dat$se.exposure,
        beta_out = dat$beta.outcome,
        beta_out_se = dat$se.outcome,
        num_gamma = 2
      )
    ))
    if (!is.null(pcmr_res$effect)) pcmr_res$effect else NA
  }, error = function(e) NA)

  res$CAUSE <- tryCatch({
    X_df <- dat %>%
      dplyr::rename(
        snp = SNP,
        beta_hat_1 = beta.exposure,
        seb1 = se.exposure,
        beta_hat_2 = beta.outcome,
        seb2 = se.outcome
      ) %>%
      dplyr::select(snp, beta_hat_1, seb1, beta_hat_2, seb2)

    if (exists("new_cause_data")) {
      X <- new_cause_data(X_df)
    } else {
      X <- X_df
      class(X) <- c("cause_data", "data.frame")
    }

    invisible(capture.output(
      suppressWarnings(suppressMessages(params <- est_cause_params(X, X$snp)))
    ))
    invisible(capture.output(
      suppressWarnings(suppressMessages(
        cause_res <- cause(X = X, variants = X$snp, param_ests = params, pval_thresh = 0.05)
      ))
    ))

    causal_quants_matrix <- cause:::summary.cause(cause_res)$quants[[2]]
    if ("gamma" %in% colnames(causal_quants_matrix)) {
      as.numeric(causal_quants_matrix[1, "gamma"])
    } else {
      NA
    }
  }, error = function(e) NA)

  res$MR_Horse <- tryCatch({
    jags_data <- list(
      by = dat$beta.outcome,
      bx = dat$beta.exposure,
      sy = dat$se.outcome,
      sx = dat$se.exposure,
      N = nrow(dat)
    )
    tmp_file <- tempfile(pattern = paste0("mr_horse_sim_", sim_id, "_"), fileext = ".txt")
    writeLines(mr_horse_model_code, tmp_file)

    fit <- R2jags::jags(
      data = jags_data,
      parameters.to.save = c("theta"),
      model.file = tmp_file,
      n.chains = 3,
      n.iter = JAGS_ITER,
      n.burnin = JAGS_BURNIN,
      progress.bar = "none"
    )
    theta_est <- fit$BUGSoutput$summary["theta", "mean"]
    if (file.exists(tmp_file)) file.remove(tmp_file)
    theta_est
  }, error = function(e) {
    if (exists("tmp_file") && file.exists(tmp_file)) file.remove(tmp_file)
    NA
  })

  res
}

# ============================================================================
# Section 5: Simulation loop
# ============================================================================
run_scenario_loop <- function(scenario_params, scenario_name, n_sim) {
  all_results <- data.frame()

  for (param in scenario_params) {
    message(sprintf("  -> Running scenario %s | severity = %g", scenario_name, param))

    sim_res <- foreach(
      i = seq_len(n_sim),
      .packages = c(
        "MASS", "dplyr", "tidyr", "MendelianRandomization",
        "MRMix", "RBMR", "MR.LDP", "mr.divw", "mr.simss",
        "mr.raps", "Matrix", "cause", "R2jags", "PCMR"
      ),
      .export = ls(envir = .GlobalEnv),
      .errorhandling = "pass"
    ) %dopar% {
      set.seed(SEED_NUMBER + as.integer(param * 1000) + i)

      data_obj <- generate_pleiotropy_data(param = param, scenario = scenario_name)
      methods_est <- apply_all_methods(data_obj$dat, data_obj$ld_matrix, sim_id = i)

      est_df <- data.frame(sim = i, param = param)
      for (m in METHOD_NAMES) {
        val <- methods_est[[m]]
        est_df[[m]] <- if (!is.null(val) && length(val) > 0) as.numeric(val[1]) else NA_real_
      }
      est_df
    }

    valid_res <- Filter(function(x) is.data.frame(x), sim_res)
    if (length(valid_res) == 0) {
      message(sprintf("[WARNING] All iterations failed for scenario %s at param = %g", scenario_name, param))
      next
    }

    raw_results <- bind_rows(valid_res)

    summary_res <- raw_results %>%
      pivot_longer(cols = any_of(METHOD_NAMES), names_to = "Method", values_to = "Estimate") %>%
      group_by(Method, param) %>%
      summarise(
        Mean_Est = mean(Estimate, na.rm = TRUE),
        SD_Est = sd(Estimate, na.rm = TRUE),
        SE_Est = ifelse(sum(!is.na(Estimate)) > 1, SD_Est / sqrt(sum(!is.na(Estimate))), 0),
        N = sum(!is.na(Estimate)),
        .groups = "drop"
      ) %>%
      filter(!is.na(Mean_Est)) %>%
      group_by(param) %>%
      mutate(Concordance_CV = sd(Mean_Est, na.rm = TRUE) / abs(mean(Mean_Est, na.rm = TRUE))) %>%
      ungroup()

    all_results <- bind_rows(all_results, summary_res)
  }

  all_results
}

# ============================================================================
# Section 6: Plotting helpers
# ============================================================================
fig2_colors <- c(
  MR_Egger   = "#000000",
  dIVW       = "#0072B2",
  MRMix      = "#D55E00",
  SimSS_IVW  = "#E69F00",
  SimSS_RAPS = "#009E73",
  MR_RAPS    = "#CC79A7",
  RBMR       = "#8C564B",
  MR_LDP     = "#9467BD",
  PCMR       = "#7F7F7F",
  CAUSE      = "#BCBD22",
  MR_Horse   = "#17BECF"
)

fig2_labels <- c(
  MR_Egger   = "MR-Egger",
  dIVW       = "dIVW",
  MRMix      = "MRMix",
  SimSS_IVW  = "MR-SimSS (IVW)",
  SimSS_RAPS = "MR-SimSS (RAPS)",
  MR_RAPS    = "MR-RAPS",
  RBMR       = "RBMR",
  MR_LDP     = "MR-LDP",
  PCMR       = "PCMR",
  CAUSE      = "CAUSE",
  MR_Horse   = "MR-Horse"
)

plot_panel_fig2 <- function(data, title, x_label) {
  ggplot(data, aes(x = param, y = Mean_Est, color = Method, group = Method)) +
    geom_hline(yintercept = TRUE_THETA, linetype = "dashed", color = "grey25", linewidth = 0.8) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    geom_errorbar(
      aes(ymin = Mean_Est - 1.96 * SE_Est, ymax = Mean_Est + 1.96 * SE_Est),
      width = 0.08,
      linewidth = 0.45
    ) +
    scale_color_manual(values = fig2_colors, labels = fig2_labels) +
    labs(
      title = title,
      x = x_label,
      y = "Estimated causal effect"
    ) +
    theme_classic(base_family = "sans", base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
}

# ============================================================================
# Section 7: Run all scenarios
# ============================================================================
message("Starting Cluster 2 simulations ...")

res_D <- run_scenario_loop(SCENARIO_PARAMS, "A", N_SIMULATIONS) # Panel D
res_E <- run_scenario_loop(SCENARIO_PARAMS, "B", N_SIMULATIONS) # Panel E
res_F <- run_scenario_loop(SCENARIO_PARAMS, "C", N_SIMULATIONS) # Panel F
res_G <- run_scenario_loop(SCENARIO_PARAMS, "D", N_SIMULATIONS) # Panel G
res_H <- run_scenario_loop(SCENARIO_PARAMS, "E", N_SIMULATIONS) # Panel H

stopCluster(cl)

method_order <- names(fig2_labels)
for (obj_name in c("res_D", "res_E", "res_F", "res_G", "res_H")) {
  tmp <- get(obj_name)
  tmp$Method <- factor(tmp$Method, levels = method_order)
  assign(obj_name, tmp)
}

# ============================================================================
# Section 8: Save summaries
# ============================================================================
write.csv(res_D, file.path(out_dir, "cluster2_panelD_summary.csv"), row.names = FALSE)
write.csv(res_E, file.path(out_dir, "cluster2_panelE_summary.csv"), row.names = FALSE)
write.csv(res_F, file.path(out_dir, "cluster2_panelF_summary.csv"), row.names = FALSE)
write.csv(res_G, file.path(out_dir, "cluster2_panelG_summary.csv"), row.names = FALSE)
write.csv(res_H, file.path(out_dir, "cluster2_panelH_summary.csv"), row.names = FALSE)

# ============================================================================
# Section 9: Assemble figure
# ============================================================================
pD <- plot_panel_fig2(res_D, "D. Outlier pleiotropy", "Outlier pleiotropy severity")
pE <- plot_panel_fig2(res_E, "E. Directional pleiotropy", "Directional pleiotropy severity")
pF <- plot_panel_fig2(res_F, "F. Balanced weak-instrument pleiotropy", "Balanced weak-instrument pleiotropy severity")
pG <- plot_panel_fig2(res_G, "G. Correlated pleiotropy (InSIDE violation)", "Correlated pleiotropy severity")
pH <- plot_panel_fig2(res_H, "H. LD-structured pleiotropy", "LD-structured pleiotropy severity")

legend_plot <- ggplot(res_D, aes(x = param, y = Mean_Est, color = Method)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = fig2_colors, labels = fig2_labels) +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

shared_legend <- get_legend(legend_plot)

panels_grid <- plot_grid(pD, pE, pF, pG, pH, NULL, ncol = 3, align = "hv")
final_plot <- plot_grid(panels_grid, shared_legend, ncol = 1, rel_heights = c(1, 0.14))

pdf_name <- "cluster2_pleiotropy_ld_figure.pdf"
png_name <- "cluster2_pleiotropy_ld_figure.png"

ggsave(
  filename = file.path(fig_dir, pdf_name),
  plot = final_plot,
  width = 17,
  height = 10,
  device = cairo_pdf
)

ggsave(
  filename = file.path(fig_dir, png_name),
  plot = final_plot,
  width = 17,
  height = 10,
  dpi = 500
)

# ============================================================================
# Section 10: Save diagnostic report
# ============================================================================
report_path <- file.path(log_dir, "cluster2_diagnostic_report.txt")
sink(report_path)

cat("======================================================\n")
cat("Cluster 2 Diagnostic Report\n")
cat("======================================================\n\n")

all_panels <- bind_rows(
  "Panel_D" = res_D,
  "Panel_E" = res_E,
  "Panel_F" = res_F,
  "Panel_G" = res_G,
  "Panel_H" = res_H,
  .id = "Panel"
)

cat("--- Mean estimates by panel, scenario parameter, and method ---\n")
print(
  all_panels %>%
    dplyr::select(Panel, param, Method, Mean_Est) %>%
    dplyr::arrange(Panel, param, Method),
  n = Inf
)

cat("\n--- Concordance (disagreement index) by panel ---\n")
print_cv <- function(df, label) {
  cat(paste0("\n", label, "\n"))
  cv_df <- df %>%
    group_by(param) %>%
    summarise(Disagreement = round(first(Concordance_CV), 4), .groups = "drop")
  print(cv_df)
}

print_cv(res_D, "Panel D: Outlier pleiotropy")
print_cv(res_E, "Panel E: Directional pleiotropy")
print_cv(res_F, "Panel F: Balanced weak-instrument pleiotropy")
print_cv(res_G, "Panel G: Correlated pleiotropy / InSIDE violation")
print_cv(res_H, "Panel H: LD-structured pleiotropy")

sink()

# ============================================================================
# Section 11: Completion message
# ============================================================================
message(
  "\n", strrep("=", 72), "\n",
  "Cluster 2 complete [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]\n",
  strrep("-", 72), "\n",
  "Summary outputs : ", out_dir, "\n",
  "Figure outputs  : ", fig_dir, "\n",
  "Log outputs     : ", log_dir, "\n",
  "Files written   :\n",
  "  ", file.path(out_dir, "cluster2_panelD_summary.csv"), "\n",
  "  ", file.path(out_dir, "cluster2_panelE_summary.csv"), "\n",
  "  ", file.path(out_dir, "cluster2_panelF_summary.csv"), "\n",
  "  ", file.path(out_dir, "cluster2_panelG_summary.csv"), "\n",
  "  ", file.path(out_dir, "cluster2_panelH_summary.csv"), "\n",
  "  ", file.path(fig_dir, pdf_name), "\n",
  "  ", file.path(fig_dir, png_name), "\n",
  "  ", report_path, "\n",
  strrep("=", 72), "\n"
)
