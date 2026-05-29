# =============================================================================
# 04b_replot_ctsem.R
# Re-plots fig_ctsem.pdf from the existing results_ctsem.rds without
# re-running any ctsem fits. Run this in the same working directory where
# results/results_ctsem.rds already exists.
# =============================================================================

source("00_helpers.R")
source("00_theme.R")

results <- readRDS("results/results_ctsem.rds")

cat("Columns in results_ctsem.rds:\n")
print(names(results))
cat(sprintf("\nNumber of replications: %d\n\n", nrow(results)))

# Quick numerical sanity check before plotting
cat("Means by estimator / method:\n")
print(
  results %>%
    summarise(
      A_ols_eta    = mean(A_ols_eta),
      A_ols_etahat = mean(A_ols_etahat),
      A_ct_eta     = mean(A_ct_eta),
      A_ct_etahat  = mean(A_ct_etahat)
    )
)


# -----------------------------------------------------------------------------
# Reshape WITHOUT regex: stack the four columns one by one and label them
# explicitly. This avoids the "eta" vs "etahat" name-collision issue.
# -----------------------------------------------------------------------------

plot_data <- bind_rows(
  data.frame(A_est = results$A_ols_eta,    estimator = "OLS (AR(1))", method = "eta"),
  data.frame(A_est = results$A_ols_etahat, estimator = "OLS (AR(1))", method = "etahat"),
  data.frame(A_est = results$A_ct_eta,     estimator = "ctsem",       method = "eta"),
  data.frame(A_est = results$A_ct_etahat,  estimator = "ctsem",       method = "etahat")
) %>%
  mutate(
    method    = factor(method,    levels = c("eta", "etahat")),
    estimator = factor(estimator, levels = c("OLS (AR(1))", "ctsem"))
  )

cat("\nRows per (estimator, method) cell:\n")
print(table(plot_data$estimator, plot_data$method))


# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------

p_ctsem <- ggplot(plot_data, aes(x = A_est, fill = method, colour = method)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  geom_vline(xintercept = -0.5,
             linetype = "dashed", linewidth = 0.4, colour = COL_REFERENCE) +
  facet_wrap(~ estimator, ncol = 1) +
  scale_fill_manual(values = THESIS_FILL, labels = THESIS_LABS) +
  scale_colour_manual(values = THESIS_FILL, labels = THESIS_LABS) +
  labs(
    x       = expression("Estimated  " * hat(A)),
    y       = "Density",
    caption = sprintf(
      "%d replications. N = 200, T = 15, dt = 1, Lambda = (0.8, 0.7, 0.9), sigma_eps = 0.5.",
      nrow(results)
    )
  ) +
  thesis_theme()

print(p_ctsem)
save_thesis_plot(p_ctsem, "results/fig_ctsem.pdf",
                 width = 4.8, height = 4.0)

cat("\nRegenerated results/fig_ctsem.pdf\n")