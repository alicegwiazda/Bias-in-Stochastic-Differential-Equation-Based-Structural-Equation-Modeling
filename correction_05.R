# =============================================================================
# 04_correction.R
# Section 8.6 (new): bias correction via estimated reliability.
#
# Joran's idea: instead of treating the bias as a phenomenon to describe, use
# the formula to remove it. The bias formula is
#
#     plim A_hat = A + log(lambda) / dt
#
# so given an estimate of lambda, we can construct
#
#     A_corrected = A_hat - log(lambda_hat) / dt
#
# and expect it to sit near the true A. The catch is that we cannot use
# the oracle lambda = cor(eta, eta_hat)^2 (it depends on the true latent).
# We need an estimate of lambda that uses only the observed indicators.
#
# Two standard choices, both single-number reliabilities of the composite:
#   - Cronbach's alpha     (psych::alpha)
#   - McDonald's omega     (psych::omega)
#
# Design: re-run the baseline (N=200, T=15, dt=1, Lambda=(0.8,0.7,0.9),
# sigma_eps=0.5), 500 replications. For each replication record:
#   - oracle lambda                 cor(eta, eta_hat)^2
#   - alpha-based lambda            psych::alpha(y1,y2,y3)$total$std.alpha
#   - omega-based lambda            psych::omega(y, 1)$omega.tot   (or omega_h)
#   - A_hat from factor scores      (already in one_replication)
#   - A_corrected with each lambda estimate
#
# Outputs:
#   results/results_correction.rds      raw per-replication results
#   results/fig_correction.pdf          three densities: biased, corrected_alpha,
#                                       corrected_omega, against true A
# =============================================================================

source("00_helpers.R")
source("00_theme.R")

set.seed(42)

dir.create("results", showWarnings = FALSE)


# -----------------------------------------------------------------------------
# Helper: estimate reliability two ways on the indicators of one replication.
# Returns a list with alpha-based lambda and omega-based lambda.
# -----------------------------------------------------------------------------
estimate_reliability <- function(data) {
  
  Y <- data[, c("y1", "y2", "y3")]
  
  # Cronbach's alpha. Standardized version, since indicators have unequal loadings.
  alpha_obj <- suppressWarnings(suppressMessages(
    psych::alpha(Y, warnings = FALSE)
  ))
  lambda_alpha <- as.numeric(alpha_obj$total$std.alpha)
  
  # McDonald's omega. We want a single number; omega.tot is the standard
  # report (total reliability of the composite under a one-factor model).
  # If you wanted the hierarchical version it would be omega_h, but tot is
  # the right analogue of "reliability of the factor score composite".
  omega_obj <- suppressWarnings(suppressMessages(
    psych::omega(Y, nfactors = 1, plot = FALSE, flip = FALSE)
  ))
  lambda_omega <- as.numeric(omega_obj$omega.tot)
  
  list(lambda_alpha = lambda_alpha, lambda_omega = lambda_omega)
}


# -----------------------------------------------------------------------------
# Per-replication function. Wraps the existing pipeline and adds the two
# reliability estimates plus the corrected A values.
# -----------------------------------------------------------------------------
one_replication_corrected <- function(
    N             = 200,
    Tpoints       = 15,
    dt            = 1,
    A_true        = -0.5,
    sigma_process = 0.5,
    loadings      = c(0.8, 0.7, 0.9),
    sigma_measure = 0.5
) {
  
  data <- simulate_ou_data(
    N             = N,
    Tpoints       = Tpoints,
    dt            = dt,
    A_true        = A_true,
    sigma_process = sigma_process,
    loadings      = loadings,
    sigma_measure = sigma_measure
  )
  data <- extract_factor_scores(data)
  
  A_from_eta    <- fit_ar1_drift(data, "eta_true", dt = dt)
  A_from_etahat <- fit_ar1_drift(data, "eta_hat",  dt = dt)
  
  if (is.na(A_from_eta) || is.na(A_from_etahat)) return(NULL)
  
  lambda_oracle <- reliability(data)
  rel           <- estimate_reliability(data)
  
  # The correction subtracts log(lambda_hat)/dt from the biased A_hat.
  # If lambda_hat == oracle lambda, the correction is exact in the limit.
  A_corr_oracle <- A_from_etahat - log(lambda_oracle)         / dt
  A_corr_alpha  <- A_from_etahat - log(rel$lambda_alpha)      / dt
  A_corr_omega  <- A_from_etahat - log(rel$lambda_omega)      / dt
  
  data.frame(
    A_true          = A_true,
    dt              = dt,
    lambda_oracle   = lambda_oracle,
    lambda_alpha    = rel$lambda_alpha,
    lambda_omega    = rel$lambda_omega,
    A_from_eta      = A_from_eta,
    A_from_etahat   = A_from_etahat,
    A_corr_oracle   = A_corr_oracle,
    A_corr_alpha    = A_corr_alpha,
    A_corr_omega    = A_corr_omega
  )
}


# -----------------------------------------------------------------------------
# Run the replications
# -----------------------------------------------------------------------------

N_REPS <- 500

cat(sprintf("Running %d replications of the correction experiment...\n", N_REPS))
results_corr <- bind_rows(lapply(seq_len(N_REPS), function(i) {
  if (i %% 50 == 0) cat(sprintf("  replication %d / %d\n", i, N_REPS))
  one_replication_corrected()
}))

saveRDS(results_corr, "results/results_correction.rds")


# -----------------------------------------------------------------------------
# Numerical summary --- the table to put in the thesis
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("CORRECTION RESULTS\n")
cat("=============================================================\n")

summary_corr <- results_corr %>%
  summarise(
    A_true             = mean(A_true),
    mean_A_etahat      = mean(A_from_etahat),
    mean_A_corr_oracle = mean(A_corr_oracle),
    mean_A_corr_alpha  = mean(A_corr_alpha),
    mean_A_corr_omega  = mean(A_corr_omega),
    bias_etahat        = mean(A_from_etahat - A_true),
    bias_corr_oracle   = mean(A_corr_oracle - A_true),
    bias_corr_alpha    = mean(A_corr_alpha  - A_true),
    bias_corr_omega    = mean(A_corr_omega  - A_true),
    rmse_etahat        = sqrt(mean((A_from_etahat - A_true)^2)),
    rmse_corr_oracle   = sqrt(mean((A_corr_oracle - A_true)^2)),
    rmse_corr_alpha    = sqrt(mean((A_corr_alpha  - A_true)^2)),
    rmse_corr_omega    = sqrt(mean((A_corr_omega  - A_true)^2)),
    mean_lambda_oracle = mean(lambda_oracle),
    mean_lambda_alpha  = mean(lambda_alpha),
    mean_lambda_omega  = mean(lambda_omega)
  )

cat(sprintf("  True A                              : %.3f\n", summary_corr$A_true))
cat("  -----------------------------------------------\n")
cat(sprintf("  Mean A from factor scores (biased)  : %.3f   bias %.3f   RMSE %.3f\n",
            summary_corr$mean_A_etahat,      summary_corr$bias_etahat,
            summary_corr$rmse_etahat))
cat(sprintf("  Corrected with oracle lambda        : %.3f   bias %.3f   RMSE %.3f\n",
            summary_corr$mean_A_corr_oracle, summary_corr$bias_corr_oracle,
            summary_corr$rmse_corr_oracle))
cat(sprintf("  Corrected with Cronbach alpha       : %.3f   bias %.3f   RMSE %.3f\n",
            summary_corr$mean_A_corr_alpha,  summary_corr$bias_corr_alpha,
            summary_corr$rmse_corr_alpha))
cat(sprintf("  Corrected with McDonald omega       : %.3f   bias %.3f   RMSE %.3f\n",
            summary_corr$mean_A_corr_omega,  summary_corr$bias_corr_omega,
            summary_corr$rmse_corr_omega))
cat("  -----------------------------------------------\n")
cat(sprintf("  Mean lambda  (oracle/alpha/omega)   : %.3f / %.3f / %.3f\n",
            summary_corr$mean_lambda_oracle,
            summary_corr$mean_lambda_alpha,
            summary_corr$mean_lambda_omega))

# -----------------------------------------------------------------------------
# Figure: density of estimates, four lines:
#   - true latents (ground-truth baseline)
#   - biased factor scores
#   - corrected with Cronbach alpha
#   - corrected with McDonald omega
# -----------------------------------------------------------------------------

plot_data <- results_corr %>%
  select(A_from_eta, A_from_etahat, A_corr_alpha, A_corr_omega) %>%
  pivot_longer(everything(),
               names_to  = "method",
               values_to = "A_est") %>%
  mutate(method = recode(method,
                         "A_from_eta"    = "true latents",
                         "A_from_etahat" = "factor scores (biased)",
                         "A_corr_alpha"  = "corrected (Cronbach alpha)",
                         "A_corr_omega"  = "corrected (McDonald omega)")) %>%
  mutate(method = factor(method, levels = c(
    "true latents",
    "factor scores (biased)",
    "corrected (Cronbach alpha)",
    "corrected (McDonald omega)"
  )))

p_corr <- ggplot(plot_data, aes(x = A_est, fill = method, colour = method)) +
  geom_density(alpha = 0.30, linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dashed",
             linewidth = 0.4, colour = COL_REFERENCE) +
  annotate("text", x = -0.5, y = 0,
           label = "true A = -0.5",
           hjust = -0.05, vjust = -0.8,
           size = 2.6, colour = "gray35", family = "serif") +
  guides(
    fill   = guide_legend(nrow = 2, byrow = TRUE),
    colour = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  labs(
    x        = expression("Estimated  " * hat(A)),
    y        = "Density",
    caption  = sprintf(
      "%d replications. N = 200, T = 15, dt = 1, loadings = (0.8, 0.7, 0.9), sigma_eps = 0.5.",
      nrow(results_corr)
    )
  ) +
  thesis_theme() +
  theme(
    legend.position  = "top",
    legend.box       = "vertical",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.4, "cm"),
    legend.spacing.y = unit(0, "pt"),
    legend.margin    = margin(t = 2, b = 2)
  )

print(p_corr)
save_thesis_plot(p_corr, "results/fig_correction.pdf",
                 width = 4.8, height = 3.2)

cat("\nFigure saved to results/fig_correction.pdf\n")