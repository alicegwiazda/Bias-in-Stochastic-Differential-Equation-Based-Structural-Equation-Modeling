# =============================================================================
# 01_baseline.R
# Section 8.1: does the predicted bias appear under default conditions?
#
# Design: N = 200, T = 15, dt = 1, A_true = -0.5,
#         loadings = (0.8, 0.7, 0.9), sigma_measure = 0.5.
# 500 replications. For each, we estimate A two ways:
#   (a) from the true latents — should recover A ~= -0.5
#   (b) from EFA factor scores — should be systematically more negative
#
# Outputs:
#   results/results_baseline.rds   raw per-replication estimates
#   results/fig_baseline.pdf       density plot of the two distributions
# =============================================================================

source("00_helpers.R")
source("00_theme.R")

set.seed(42)

dir.create("results", showWarnings = FALSE)


# -----------------------------------------------------------------------------
# Run the replications
# -----------------------------------------------------------------------------

N_REPS <- 500

cat(sprintf("Running %d replications of the baseline simulation...\n", N_REPS))
results <- bind_rows(lapply(seq_len(N_REPS), function(i) {
  if (i %% 50 == 0) cat(sprintf("  replication %d / %d\n", i, N_REPS))
  one_replication()
}))

saveRDS(results, "results/results_baseline.rds")
cat(sprintf("\nSaved %d successful replications to results/results_baseline.rds\n",
            nrow(results)))


# -----------------------------------------------------------------------------
# Numerical summary — copy these numbers into section 8.1 of the thesis
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("BASELINE RESULTS  (section 8.1)\n")
cat("=============================================================\n")

summary_table <- results %>%
  summarise(
    mean_A_eta     = mean(A_from_eta),
    mean_A_etahat  = mean(A_from_etahat),
    mean_bias_eta  = mean(bias_eta),
    mean_bias_hat  = mean(bias_etahat),
    sd_A_eta       = sd(A_from_eta),
    sd_A_etahat    = sd(A_from_etahat),
    rmse_eta       = sqrt(mean(bias_eta^2)),
    rmse_etahat    = sqrt(mean(bias_etahat^2)),
    predicted_bias = mean(predicted_bias)
  )

cat(sprintf("  True A                            : %.3f\n", results$A_true[1]))
cat(sprintf("  Mean A from true latents          : %.3f  (sd %.3f)\n",
            summary_table$mean_A_eta, summary_table$sd_A_eta))
cat(sprintf("  Mean A from EFA factor scores     : %.3f  (sd %.3f)\n",
            summary_table$mean_A_etahat, summary_table$sd_A_etahat))
cat(sprintf("  Empirical bias (factor scores)    : %.3f\n",
            summary_table$mean_bias_hat))
cat(sprintf("  Analytical prediction log(lambda)/dt: %.3f\n",
            summary_table$predicted_bias))
cat(sprintf("  RMSE (true latents)               : %.3f\n",
            summary_table$rmse_eta))
cat(sprintf("  RMSE (factor scores)              : %.3f\n",
            summary_table$rmse_etahat))


# -----------------------------------------------------------------------------
# Figure: density of the two distributions
# -----------------------------------------------------------------------------

plot_data <- results %>%
  select(A_from_eta, A_from_etahat) %>%
  pivot_longer(everything(),
               names_to  = "method",
               values_to = "A_est") %>%
  mutate(method = recode(method,
                         "A_from_eta"    = "eta",
                         "A_from_etahat" = "etahat"))

p_baseline <- ggplot(plot_data, aes(x = A_est, fill = method, colour = method)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  geom_vline(xintercept = -0.5,
             linetype = "dashed", linewidth = 0.4, colour = COL_REFERENCE) +
  annotate("text", x = -0.5, y = 0,
           label = "true A = -0.5",
           hjust = -0.05, vjust = -0.8,
           size = 2.6, colour = "gray35", family = "serif") +
  scale_fill_manual(values = THESIS_FILL, labels = THESIS_LABS) +
  scale_colour_manual(values = THESIS_FILL, labels = THESIS_LABS) +
  labs(
    x        = expression("Estimated  " * hat(A)),
    y        = "Density",
    caption  = sprintf(
      "%d replications. N = 200, T = 15, dt = 1, loadings = (0.8, 0.7, 0.9), sigma_eps = 0.5.",
      nrow(results)
    )
  ) +
  thesis_theme()

print(p_baseline)

save_thesis_plot(p_baseline, "results/fig_baseline.pdf",
                 width = 3.4, height = 2.4)

cat("\nFigure saved to results/fig_baseline.pdf\n")