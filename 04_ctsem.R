# =============================================================================
# 04_ctsem.R
# Section 8.4: validation in ctsem.
#
# Sections 8.1-8.3 estimated A via OLS on the AR(1) recursion, then
# A_hat = log(B_hat) / dt. The reduction argument in Section 7.1 claims
# that under the two-step configuration (Lambda = [1], MANIFESTVAR = 0,
# equal intervals), CTSEM computes the same estimand. This section
# confirms that empirically: we re-run the baseline (Section 8.1) using
# the ctsem package and check that the bias matches.
#
# Each replication fits two ctsem models: one on the true latent
# trajectory, one on the EFA factor-score trajectory. ctsem is much
# slower than OLS, so we use fewer replications (100) and a single
# baseline condition rather than a sweep.
#
# Outputs:
#   results/results_ctsem.rds   raw per-replication estimates
#   results/fig_ctsem.pdf       density of estimates, AR(1) vs ctsem
# =============================================================================

source("00_helpers.R")
source("00_theme.R")

suppressPackageStartupMessages({
  library(ctsem)
})

set.seed(42)

dir.create("results", showWarnings = FALSE)


# -----------------------------------------------------------------------------
# fit_ctsem_one_var()
#   Fits a one-latent, one-manifest CTSEM with LAMBDA = [1] and manifest
#   variance fixed at near-zero. Under this configuration, ctsem treats the
#   input column as the latent state directly and the Kalman filter passes
#   it through unchanged, so the MLE is equivalent to OLS on the AR(1)
#   recursion (see Section 7.1). Returns the estimated drift A, or NA on
#   failure.
# -----------------------------------------------------------------------------

fit_ctsem_one_var <- function(data, varname) {
  
  ctdata <- data %>%
    dplyr::select(id, time, dplyr::all_of(varname))
  names(ctdata)[3] <- "Y"
  ctdata$id   <- as.numeric(as.factor(ctdata$id))
  ctdata$time <- as.numeric(ctdata$time)
  
  model <- ctModel(
    type          = "ct",
    n.latent      = 1,
    n.manifest    = 1,
    manifestNames = "Y",
    latentNames   = "eta",
    LAMBDA        = matrix(1, 1, 1),
    DRIFT         = matrix("drift",     1, 1),
    DIFFUSION     = matrix("diffusion", 1, 1),
    CINT          = matrix(0, 1, 1),
    T0MEANS       = matrix(0, 1, 1),
    T0VAR         = matrix("t0var", 1, 1),
    MANIFESTMEANS = matrix(0,    1, 1),
    MANIFESTVAR   = matrix(1e-6, 1, 1)
  )
  
  fit <- tryCatch({
    suppressMessages(suppressWarnings(
      ctStanFit(
        datalong    = ctdata,
        ctstanmodel = model,
        optimize    = TRUE,
        priors      = FALSE,
        verbose     = 0
      )
    ))
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(NA_real_)
  
  as.numeric(summary(fit)$popmeans["drift", "mean"])
}


# -----------------------------------------------------------------------------
# one_ctsem_replication()
#   Same pipeline as one_replication() in 00_helpers, but uses ctsem instead
#   of OLS for the dynamic step. Returns a row with both the ctsem estimates
#   and the OLS estimates so we can compare them directly.
# -----------------------------------------------------------------------------

one_ctsem_replication <- function(
    N             = 200,
    Tpoints       = 15,
    dt            = 1,
    A_true        = -0.5,
    sigma_process = 0.5,
    loadings      = c(0.8, 0.7, 0.9),
    sigma_measure = 0.5
) {
  
  data <- simulate_ou_data(
    N             = N, Tpoints = Tpoints, dt = dt,
    A_true        = A_true, sigma_process = sigma_process,
    loadings      = loadings, sigma_measure = sigma_measure
  )
  data <- extract_factor_scores(data)
  
  # OLS estimates (the AR(1) path used in 8.1-8.3)
  A_ols_eta    <- fit_ar1_drift(data, "eta_true", dt = dt)
  A_ols_etahat <- fit_ar1_drift(data, "eta_hat",  dt = dt)
  
  # ctsem estimates (the new path)
  A_ct_eta    <- fit_ctsem_one_var(data, "eta_true")
  A_ct_etahat <- fit_ctsem_one_var(data, "eta_hat")
  
  if (anyNA(c(A_ols_eta, A_ols_etahat, A_ct_eta, A_ct_etahat))) return(NULL)
  
  data.frame(
    A_true       = A_true,
    lambda       = reliability(data),
    A_ols_eta    = A_ols_eta,
    A_ols_etahat = A_ols_etahat,
    A_ct_eta     = A_ct_eta,
    A_ct_etahat  = A_ct_etahat,
    bias_ols_eta    = A_ols_eta    - A_true,
    bias_ols_etahat = A_ols_etahat - A_true,
    bias_ct_eta     = A_ct_eta     - A_true,
    bias_ct_etahat  = A_ct_etahat  - A_true
  )
}


# -----------------------------------------------------------------------------
# Run the replications
# -----------------------------------------------------------------------------

N_REPS <- 100   # ctsem is slow; expect roughly 30-90 minutes for 100 reps

cat(sprintf("Running %d ctsem replications (this is slow)...\n", N_REPS))
results <- bind_rows(lapply(seq_len(N_REPS), function(i) {
  if (i %% 5 == 0) cat(sprintf("  replication %d / %d\n", i, N_REPS))
  one_ctsem_replication()
}))

saveRDS(results, "results/results_ctsem.rds")
cat(sprintf("\nSaved %d successful replications.\n", nrow(results)))


# -----------------------------------------------------------------------------
# Numerical summary
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("CTSEM VALIDATION RESULTS\n")
cat("=============================================================\n")
summary_ct <- results %>%
  summarise(
    mean_A_ols_eta      = mean(A_ols_eta),
    mean_A_ols_etahat   = mean(A_ols_etahat),
    mean_A_ct_eta       = mean(A_ct_eta),
    mean_A_ct_etahat    = mean(A_ct_etahat),
    sd_A_ols_etahat     = sd(A_ols_etahat),
    sd_A_ct_etahat      = sd(A_ct_etahat),
    mean_bias_ols_etahat = mean(bias_ols_etahat),
    mean_bias_ct_etahat  = mean(bias_ct_etahat),
    mean_lambda          = mean(lambda)
  )

cat(sprintf("  True A                                : %.3f\n", results$A_true[1]))
cat(sprintf("  Mean A (OLS, true latents)            : %.3f\n",
            summary_ct$mean_A_ols_eta))
cat(sprintf("  Mean A (ctsem, true latents)          : %.3f\n",
            summary_ct$mean_A_ct_eta))
cat(sprintf("  Mean A (OLS, factor scores)           : %.3f  (sd %.3f)\n",
            summary_ct$mean_A_ols_etahat, summary_ct$sd_A_ols_etahat))
cat(sprintf("  Mean A (ctsem, factor scores)         : %.3f  (sd %.3f)\n",
            summary_ct$mean_A_ct_etahat, summary_ct$sd_A_ct_etahat))
cat(sprintf("  Bias (OLS, factor scores)             : %.3f\n",
            summary_ct$mean_bias_ols_etahat))
cat(sprintf("  Bias (ctsem, factor scores)           : %.3f\n",
            summary_ct$mean_bias_ct_etahat))
cat(sprintf("  Mean lambda                           : %.3f\n",
            summary_ct$mean_lambda))


# -----------------------------------------------------------------------------
# Figure: densities of the four estimate distributions
# -----------------------------------------------------------------------------

plot_data <- results %>%
  select(A_ols_eta, A_ols_etahat, A_ct_eta, A_ct_etahat) %>%
  pivot_longer(everything(),
               names_to       = c("estimator", "method"),
               names_pattern  = "A_(ols|ct)_(eta|etahat)",
               values_to      = "A_est") %>%
  mutate(
    method    = factor(method,    levels = c("eta", "etahat")),
    estimator = factor(estimator,
                       levels = c("ols", "ct"),
                       labels = c("OLS (AR(1))", "ctsem"))
  )

p_ctsem <- ggplot(plot_data, aes(x = A_est, fill = method, colour = method)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  geom_vline(xintercept = -0.5,
             linetype = "dashed", linewidth = 0.4, colour = COL_REFERENCE) +
  facet_wrap(~ estimator, ncol = 1) +
  scale_fill_manual(values = THESIS_FILL, labels = THESIS_LABS) +
  scale_colour_manual(values = THESIS_FILL, labels = THESIS_LABS) +
  labs(
    x        = expression("Estimated  " * hat(A)),
    y        = "Density",
    caption  = sprintf(
      "%d replications. N = 200, T = 15, dt = 1, Lambda = (0.8, 0.7, 0.9), sigma_eps = 0.5.",
      nrow(results)
    )
  ) +
  thesis_theme()

print(p_ctsem)
save_thesis_plot(p_ctsem, "results/fig_ctsem.pdf",
                 width = 4.8, height = 4.0)

cat("\nFigure saved to results/fig_ctsem.pdf\n")