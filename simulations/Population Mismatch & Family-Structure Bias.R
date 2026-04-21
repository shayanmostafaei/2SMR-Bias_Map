###############################################################################
# Bias Cluster 3: Population Mismatch & Family-Structure Bias
###############################################################################
rm(list = ls())

# 0. GLOBAL CONFIGURATION
N_SIMULATIONS   <- 50       
SCENARIO_PARAMS <- seq(0, 1, by = 0.25) 
N_SNPS          <- 500 
TRUE_THETA      <- 0.3       
SAMPLE_SIZE_N1  <- 50000     
SAMPLE_SIZE_N2  <- 50000     
BASE_VARIANCE   <- 0.05      

NUM_CORES       <- 4        
SEED_NUMBER     <- 20260414  

# --- MR-Horse Specific Settings ---
JAGS_ITER       <- 10000    # Iterations for JAGS model
JAGS_BURNIN     <- 5000     # Burn-in for JAGS model   

suppressPackageStartupMessages({
  library(MASS)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(MendelianRandomization)
  library(mr.raps)
  library(mr.simss)
  library(doParallel)
  library(foreach)
  library(meta)           
  library(MRSamePopTest) 
  library(MRMix)
  library(PCMR)
  library(R2jags)
})

set.seed(SEED_NUMBER)

if (exists("cl")) { try(stopCluster(cl), silent = TRUE) }
cl <- makeCluster(NUM_CORES)
registerDoParallel(cl)

###############################################################################
# 1. DATA GENERATION MODEL & LD MATRIX
###############################################################################

# Function to create LD matrix based on Toeplitz blocks
create_ld_matrix <- function(n_snps, block_size = 10, rho = 0.5) {
  n_blocks <- ceiling(n_snps / block_size)
  base_block <- toeplitz(rho^(0:(block_size - 1)))
  ld_mat <- as.matrix(Matrix::bdiag(replicate(n_blocks, base_block, simplify = FALSE)))
  return(ld_mat[1:n_snps, 1:n_snps])
}

generate_cluster3_data <- function(param, scenario = "Population_Mismatch", n_snps = N_SNPS) {
  
  # Common base parameters
   maf_exp <- rbeta(n_snps, 1, 3) * 0.4 + 0.05
   
   drift_sd <- param * 0.1 
   maf_drift <- rnorm(n_snps, mean = 0, sd = drift_sd)

   maf_out <- maf_exp + maf_drift

   maf_out <- pmax(0.01, pmin(0.99, maf_out))

# ---------------------------------------------------------
se_exp <- sqrt(1 / (2 * maf_exp * (1 - maf_exp) * SAMPLE_SIZE_N1))
se_out <- sqrt(1 / (2 * maf_out * (1 - maf_out) * SAMPLE_SIZE_N2))


  
  beta_exp_pop2 <- rep(NA, n_snps) # Initial NA value for second population
  
  if (scenario == "Population_Mismatch") {
    # Panel A: Two different LD matrices
    rho_exp <- 0.5
    rho_out <- max(0.1, 0.5 - (param * 0.4)) 
    
    ld_matrix_exp <- create_ld_matrix(n_snps, rho = rho_exp)
    ld_matrix_out <- create_ld_matrix(n_snps, rho = rho_out)
    
    alpha_exp <- MASS::mvrnorm(1, mu = rep(0, n_snps), Sigma = ld_matrix_exp * BASE_VARIANCE)
    alpha_out <- MASS::mvrnorm(1, mu = rep(0, n_snps), Sigma = ld_matrix_out * BASE_VARIANCE)
    
    shrinkage_factor <- 1 - (param * 0.5)
    mismatch_noise <- rnorm(n_snps, mean = param * 0.02, sd = param * sqrt(BASE_VARIANCE) * 0.5)
    alpha_pop2 <- (alpha_out * shrinkage_factor) + mismatch_noise + (param * 0.05)
    
    beta_exp <- alpha_exp + rnorm(n_snps, 0, se_exp)
    beta_out <- (alpha_pop2 * TRUE_THETA) + rnorm(n_snps, 0, se_out)
    beta_exp_pop2 <- alpha_pop2 + rnorm(n_snps, 0, se_exp)
    
    se_exp_wf <- se_exp * 1.5
    se_out_wf <- se_out * 1.5
    beta_exp_wf <- alpha_exp + rnorm(n_snps, 0, se_exp_wf)
    beta_out_wf <- (alpha_pop2 * TRUE_THETA) + rnorm(n_snps, 0, se_out_wf)
    
    ld_matrix_to_return <- ld_matrix_exp
    
  } else if (scenario == "Family_Structure") {
    # Panel B: A shared LD matrix
    rho_shared <- 0.5
    ld_matrix_shared <- create_ld_matrix(n_snps, rho = rho_shared)
    
    alpha <- MASS::mvrnorm(1, mu = rep(0, n_snps), Sigma = ld_matrix_shared * BASE_VARIANCE)
    
    se_exp_wf <- se_exp * (1 + param * 2.0) 
    se_out_wf <- se_out * (1 + param * 2.0)
    
    dynastic_confounding <- (alpha * param * 0.8) + rnorm(n_snps, 0.02, param * 0.05)
    
    beta_exp <- alpha + dynastic_confounding + rnorm(n_snps, 0, se_exp)
    beta_out <- (alpha * TRUE_THETA) + dynastic_confounding + rnorm(n_snps, 0, se_out)
    
    beta_exp_wf <- alpha + rnorm(n_snps, 0, se_exp_wf)
    beta_out_wf <- (alpha * TRUE_THETA) + rnorm(n_snps, 0, se_out_wf)
    
    ld_matrix_to_return <- ld_matrix_shared
  }
  
  dat <- data.frame(
    SNP = paste0("rs", 1:n_snps),
    beta.exposure = beta_exp,
    se.exposure = se_exp,
    beta.outcome = beta_out,
    se.outcome = se_out,
    beta.exp_pop2 = beta_exp_pop2, 
    se.exp_pop2 = se_out,          
    beta.exposure.wf = beta_exp_wf,
    se.exposure.wf = se_exp_wf,
    beta.outcome.wf = beta_out_wf,
    se.outcome.wf = se_out_wf
  )
  
  return(list(dat = dat, ld_matrix = ld_matrix_to_return))
}

###############################################################################
# 2. HELPER MODELS & PATCHES
###############################################################################

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

METHOD_NAMES <- c("MR_Egger", "MR_RAPS", 
                  "IVW_SamePopFilter", "IVW_WithinFamily", 
                  "SimSS_IVW", "SimSS_RAPS",
                  "MRMix", "PCMR", "MR_Horse")

###############################################################################
# 3. ESTIMATION WRAPPER
###############################################################################

apply_cluster3_methods <- function(data_list, sim_id) {
  dat <- data_list$dat
  ld_matrix <- data_list$ld_matrix
  
  res <- list()
  mr_obj <- tryCatch(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                              by = dat$beta.outcome, byse = dat$se.outcome), 
                     error = function(e) NULL)
  if (is.null(mr_obj)) return(res)
  
  # Classical methods
  res$MR_Egger <- tryCatch(mr_egger(mr_obj)$Estimate, error=function(e) NA)
  res$MR_RAPS <- tryCatch(mr.raps(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, over.dispersion = TRUE)$beta.hat, error=function(e) NA)
  
  res$IVW_SamePopFilter <- tryCatch({
    if (all(is.na(dat$beta.exp_pop2))) { NA_real_ } else {
      pdf(NULL); spt_res <- same_pop_test(Bexp = dat$beta.exposure, Bout = dat$beta.exp_pop2, SEexp = dat$se.exposure, SEout = dat$se.exp_pop2, SNPlist = dat$SNP, Fisher = FALSE); dev.off()
      valid_idx <- spt_res$pval > 0.05 
      if (sum(valid_idx) > 3) {
        mr_ivw(mr_input(bx = dat$beta.exposure[valid_idx], bxse = dat$se.exposure[valid_idx], by = dat$beta.outcome[valid_idx], byse = dat$se.outcome[valid_idx]))$Estimate
      } else { NA_real_ }
    }
  }, error=function(e) NA)
  
  res$IVW_WithinFamily <- tryCatch({
    mr_ivw(mr_input(bx = dat$beta.exposure.wf, bxse = dat$se.exposure.wf, by = dat$beta.outcome.wf, byse = dat$se.outcome.wf))$Estimate
  }, error=function(e) NA)
  
  pval_thresh <- pchisq(10, df=1, lower.tail=FALSE) 
  res$SimSS_IVW <- tryCatch(mr_simss(dat, mr_method="mr_ivw", threshold=pval_thresh, splits=2, n.iter=100)$summary$b, error=function(e) NA)
  res$SimSS_RAPS <- tryCatch(mr_simss(dat, mr_method="mr_raps", threshold=pval_thresh, splits=2, n.iter=100)$summary$b, error=function(e) NA)
  
  # --- Advanced / New models ---
  res$MRMix <- tryCatch(MRMix(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)$theta, error=function(e) NA)
  
  res$PCMR <- tryCatch({
    suppressWarnings(capture.output(pcmr_res <- PCMR(beta_ex = dat$beta.exposure, beta_ex_se = dat$se.exposure,
                                                    beta_out = dat$beta.outcome, beta_out_se = dat$se.outcome, num_gamma = 2)))
    if (!is.null(pcmr_res$effect)) pcmr_res$effect else NA
  }, error=function(e) NA)

  res$MR_Horse <- tryCatch({
    jags_data <- list(by = dat$beta.outcome, bx = dat$beta.exposure, sy = dat$se.outcome, sx = dat$se.exposure, N = nrow(dat))
    tmp_file <- tempfile(pattern = paste0("mr_horse_sim_", sim_id, "_"), fileext = ".txt")
    writeLines(mr_horse_model_code, tmp_file)
    fit <- R2jags::jags(data = jags_data, parameters.to.save = c("theta"), model.file = tmp_file,
                        n.chains = 3, n.iter = JAGS_ITER, n.burnin = JAGS_BURNIN, progress.bar = "none")
    theta_est <- fit$BUGSoutput$summary["theta", "mean"]
    if (file.exists(tmp_file)) file.remove(tmp_file)
    theta_est
  }, error = function(e) { if (exists("tmp_file") && file.exists(tmp_file)) file.remove(tmp_file); NA })
  
  return(res)
}

###############################################################################
# 4. SIMULATION LOOP (PARALLELIZED, WARNINGS, CV & FILE SAVING)
###############################################################################

# Initialize a global list to capture warnings across all scenarios
GLOBAL_WARNING_LOG <- list()

run_scenario_loop <- function(scenario_params, scenario_name, n_sim) {
  all_results <- data.frame()
  
  for (param in scenario_params) {
    cat(sprintf("  -> Running %s | param = %g\n", scenario_name, param))
    
    sim_res <- foreach(i = 1:n_sim, 
                       .packages = c("MASS", "Matrix", "dplyr", "tidyr", "MendelianRandomization", 
                                     "mr.raps", "meta", "MRSamePopTest", "mr.simss",
                                     "MRMix", "PCMR", "R2jags"),
                       .export = ls(envir = .GlobalEnv),
                       .errorhandling = "pass") %dopar% {
                         
                         set.seed(SEED_NUMBER + param * 1000 + i) 
                         data_list <- generate_cluster3_data(param = param, scenario = scenario_name)
                         methods_est <- apply_cluster3_methods(data_list, sim_id = i)
                         
                         est_df <- data.frame(sim = i, param = param)
                         for (m in METHOD_NAMES) {
                           val <- methods_est[[m]]
                           est_df[[m]] <- if (!is.null(val) && length(val) > 0) as.numeric(val[1]) else NA_real_
                         }
                         est_df
                       }
    
    valid_res <- Filter(function(x) is.data.frame(x), sim_res)
    
    summary_res <- bind_rows(valid_res) %>%
      pivot_longer(cols = any_of(METHOD_NAMES), names_to = "Method", values_to = "Estimate") %>%
      group_by(Method, param) %>%
      summarise(
        Mean_Est = mean(Estimate, na.rm = TRUE),
        SE_Est = sd(Estimate, na.rm = TRUE) / sqrt(max(sum(!is.na(Estimate)), 1)),
        .groups = 'drop'
      ) %>%
      group_by(param) %>%
      mutate(Concordance_CV = sd(Mean_Est, na.rm = TRUE) / abs(mean(Mean_Est, na.rm = TRUE))) %>%
      ungroup()
    
    # Log and print warnings for models that returned NA
    na_models <- summary_res %>% filter(is.na(Mean_Est) | is.nan(Mean_Est))
    if(nrow(na_models) > 0) {
      for(m in unique(na_models$Method)) {
        warn_msg <- sprintf("[WARNING] Method '%s' failed (returned NA) in scenario '%s' at param = %g.", m, scenario_name, param)
        cat("     ", warn_msg, "\n")
        # Save warning to the global log using <<-
        GLOBAL_WARNING_LOG <<- append(GLOBAL_WARNING_LOG, warn_msg)
      }
    }
    
    all_results <- rbind(all_results, summary_res)
  }
  
  return(all_results)
}

###############################################################################
# 5. PLOTTING & EXECUTION
###############################################################################

my_colors <- c(                 
  "MR_Egger" = "#E31A1C",             
  "MR_RAPS" = "#CC79A7",              
  "IVW_SamePopFilter" = "#8B4513",    
  "IVW_WithinFamily" = "#800080",     
  "SimSS_IVW" = "#E69F00",            
  "SimSS_RAPS" = "#0072B2",
  "MRMix" = "#FF00FF",     
  "PCMR" = "#555555",      
  "MR_Horse" = "#000000"   
)

my_labels <- c(
  "MR_Egger" = "1. MR-Egger", 
  "MR_RAPS" = "2. MR-RAPS", 
  "IVW_SamePopFilter" = "3. IVW (Pop Filtered)", 
  "IVW_WithinFamily" = "4. Within-Family IVW",
  "SimSS_IVW" = "5. MR-SimSS (IVW)", 
  "SimSS_RAPS" = "6. MR-SimSS (RAPS)",
  "MRMix" = "7. MRMix",
  "PCMR" = "8. PCMR",
  "MR_Horse" = "9. MR-Horse"
)

plot_cluster3 <- function(data, title, x_lab) {
  data$Method <- factor(data$Method, levels = METHOD_NAMES)
  
  ggplot(data, aes(x = param, y = Mean_Est, color = Method, group = Method)) +
    geom_hline(yintercept = TRUE_THETA, linetype = "dashed", color = "darkgreen", linewidth = 1.2) +
    geom_line(linewidth = 1.2, na.rm = TRUE) + 
    geom_point(size = 3, na.rm = TRUE) +
    geom_errorbar(aes(ymin = Mean_Est - 1.96*SE_Est, ymax = Mean_Est + 1.96*SE_Est), width = 0.03, linewidth = 0.8, na.rm = TRUE) +
    # 2. Add 'breaks = METHOD_NAMES' to strictly enforce the legend order
    scale_color_manual(values = my_colors, labels = my_labels, breaks = METHOD_NAMES) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 14)) +
    labs(title = title, x = x_lab, y = "Estimated Causal Effect")
}


cat(sprintf("\nStarting Bias Cluster 3 Simulations (%d reps)...\n", N_SIMULATIONS))

res_A <- run_scenario_loop(SCENARIO_PARAMS, "Population_Mismatch", N_SIMULATIONS)
cat("\n---------------------------------------------------\n")
res_B <- run_scenario_loop(SCENARIO_PARAMS, "Family_Structure", N_SIMULATIONS)

stopCluster(cl)

pA <- plot_cluster3(res_A, "A) Population Mismatch Bias (Design Violation)", "Degree of Population Divergence (param)")
pB <- plot_cluster3(res_B, "B) Family-Structure Bias (Dynastic Effects)", "Strength of Dynastic Confounding (param)")

shared_legend <- get_legend(pA + guides(color = guide_legend(nrow = 2, byrow = TRUE)))
panels_grid <- plot_grid(
  pA + theme(legend.position="none"), 
  pB + theme(legend.position="none"), 
  ncol = 2
)
final_plot <- plot_grid(panels_grid, shared_legend, ncol = 1, rel_heights = c(1, 0.15))

file_name = "Fig2_Bias_Cluster3_Test9"
pdf_name <- paste0(file_name, ".pdf")
png_name <- paste0(file_name, ".png")

suppressWarnings({
  ggsave(pdf_name, final_plot, width = 17, height = 11, device = cairo_pdf)
  ggsave(png_name, final_plot, width = 17, height = 11, dpi = 300)
})
cat(sprintf("\nSimulation Complete! Plots saved as %s and %s\n", pdf_name, png_name))

# ====================================================================
# 6. SAVE DIAGNOSTICS AND REPORTS TO A TEXT FILE
# ====================================================================
txt_filename <- paste0(file_name, ".txt")

sink(txt_filename)

cat("======================================================\n")
cat("                DIAGNOSTIC REPORT (CLUSTER 3)\n")
cat("======================================================\n\n")

cat("--- Method Counts in res_A ---\n")
print(table(res_A$Method))
cat("\n")

cat("--- Method Counts in res_B ---\n")
print(table(res_B$Method))
cat("\n")

cat("--- Missing Values in res_A ---\n")
missing_vals_A <- res_A %>%
  group_by(Method) %>%
  summarise(Missing_Values = sum(is.na(Mean_Est)))
print(as.data.frame(missing_vals_A))
cat("\n")

cat("--- Missing Values in res_B ---\n")
missing_vals_B <- res_B %>%
  group_by(Method) %>%
  summarise(Missing_Values = sum(is.na(Mean_Est)))
print(as.data.frame(missing_vals_B))
cat("\n")

cat("======================================================\n")
cat("                ALL PANELS DATA (Mean Estimates)\n")
cat("======================================================\n\n")

all_panels <- dplyr::bind_rows(
  "Panel_A" = res_A,
  "Panel_B" = res_B,
  .id = "Panel"
)

all_panels_summary <- all_panels %>%
  dplyr::select(Panel, param, Method, Mean_Est) %>%
  dplyr::arrange(Panel, param, Method)

print(all_panels_summary, n = Inf)

cat("\n======================================================\n")
cat(" FINAL CONCORDANCE REPORT (Disagreement Index) \n")
cat("======================================================\n")

print_cv <- function(df, label) {
  cat(paste("\nScenario", label, ":\n"))
  cv_df <- df %>% 
    group_by(param) %>% 
    summarise(Disagreement = round(first(Concordance_CV), 4), .groups = "drop")
  print(as.data.frame(cv_df))
}

print_cv(res_A, "A (Population Mismatch)")
print_cv(res_B, "B (Family Structure)")

# ---------------------------------------------------------
# APPEND WARNINGS AND ERRORS
# ---------------------------------------------------------
cat("\n======================================================\n")
cat("               WARNINGS & ERROR LOGS \n")
cat("======================================================\n")

if (length(GLOBAL_WARNING_LOG) > 0) {
  for (warn in GLOBAL_WARNING_LOG) {
    cat(warn, "\n")
  }
} else {
  cat("No errors or warnings recorded during the simulation.\n")
}

sink()

cat(sprintf("\n[SUCCESS] Diagnostic, Concordance, and Warning logs have been successfully saved to: %s/%s\n", getwd(), txt_filename))
