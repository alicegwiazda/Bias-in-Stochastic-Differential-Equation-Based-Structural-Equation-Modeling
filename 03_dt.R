# =============================================================================
# 03_dt.R
# Section 8.3: bias as a function of the observation interval.
#
# This is the CT-specific prediction. The bias formula log(lambda)/dt scales
# inversely with dt: halving dt should roughly double the (negative) bias,
# holding lambda fixed. To test this we hold the measurement model fixed
# (Lambda = (0.8, 0.7, 0.9), sigma_eps = 0.5) and vary dt only.
#
# A subtlety: the AR coefficient B = exp(A * dt) also depends on dt. To keep
# the underlying continuous-time process the same across conditions, we hold
# A_true and sigma_process fixed and let B vary with dt. This is the right
# thing to do --- we are simulating one continuous process sampled at different
# rates, not multiple processes.
#
# Outputs:
#   results/results_dt.rds      raw per-replication results
#   results/fig_dt.pdf          empirical bias vs analytical prediction
# =============================================================================

source("00_helpers.R")
source("00_theme.R")

set.seed(42)

dir.create("results", showWarnings = FALSE)


# -----------------------------------------------------------------------------
# Design
# -----------------------------------------------------------------------------

N_REPS_PER_COND <- 100

# Vary dt while holding everything else fixed.
# We pick values that span an order of magnitude to make the 1/dt
# scaling visually unambiguous: a small dt (frequent sampling) should
# give a much larger bias than a large dt (sparse sampling).
DT_GRID <- c(0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0)


# -----------------------------------------------------------------------------
# Run the replications
# -----------------------------------------------------------------------------

cat("Running dt study (100 reps per condition)...\n")
results_dt <- bind_rows(lapply(DT_GRID, function(dt_val) {
  cat(sprintf("  dt = %.2f\n", dt_val))
  bind_rows(lapply(seq_len(N_REPS_PER_COND), function(i) {
    one_replication(dt = dt_val)
  }))
}))

saveRDS(results_dt, "results/results_dt.rds")


# -----------------------------------------------------------------------------
# Numerical summary
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("BIAS BY OBSERVATION INTERVAL\n")
cat("=============================================================\n")
print(
  results_dt %>%
    group_by(dt) %>%
    summarise(
      lambda             = round(mean(reliability),    3),
      bias_true_latents  = round(mean(bias_eta),       3),
      bias_factor_scores = round(mean(bias_etahat),    3),
      predicted_bias     = round(mean(predicted_bias), 3),
      n                  = n()
    )
)


# -----------------------------------------------------------------------------
# Figure: empirical bias vs analytical prediction, both as a function of dt
# -----------------------------------------------------------------------------

agg <- results_dt %>%
  group_by(dt) %>%
  summarise(
    lambda             = mean(reliability),
    mean_bias_eta      = mean(bias_eta),
    se_bias_eta        = sd(bias_eta) / sqrt(n()),
    mean_bias_etahat   = mean(bias_etahat),
    se_bias_etahat     = sd(bias_etahat) / sqrt(n()),
    predicted_bias     = mean(predicted_bias),
    .groups            = "drop"
  ) %>%
  pivot_longer(
    cols          = starts_with("mean_bias") | starts_with("se_bias"),
    names_to      = c(".value", "method"),
    names_pattern = "(mean_bias|se_bias)_(.*)"
  )

# Build a smooth analytical curve over a fine grid of dt values, using the
# mean lambda achieved across all conditions. Because the measurement model
# is identical across dt conditions, lambda is roughly (but not exactly)
# constant; we use its average value across the grid for a single honest
# prediction line. Any vertical gap between the empirical points and this
# line therefore reflects either (a) small variation in lambda across
# conditions, or (b) higher-order effects not captured by the leading
# log(lambda)/dt term.

mean_lambda <- mean(results_dt$reliability)
dt_fine     <- seq(min(DT_GRID), max(DT_GRID), length.out = 200)
prediction_curve <- data.frame(
  dt             = dt_fine,
  predicted_bias = log(mean_lambda) / dt_fine
)

p_dt <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted",
             linewidth = 0.4, colour = COL_REFERENCE) +
  geom_line(
    data = prediction_curve,
    aes(x = dt, y = predicted_bias, linetype = "predicted"),
    colour = COL_PREDICTED, linewidth = 0.55
  ) +
  geom_errorbar(
    data = agg,
    aes(x = dt, ymin = mean_bias - 1.96 * se_bias,
        ymax = mean_bias + 1.96 * se_bias,
        colour = method),
    width = 0.06, linewidth = 0.35
  ) +
  geom_line(
    data = agg,
    aes(x = dt, y = mean_bias, colour = method),
    linewidth = 0.5
  ) +
  geom_point(
    data = agg,
    aes(x = dt, y = mean_bias, colour = method),
    size = 1.9, shape = 16
  ) +
  scale_colour_manual(
    name   = NULL,
    values = c("eta" = COL_ETA, "etahat" = COL_ETAHAT),
    labels = THESIS_LABS
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c("predicted" = "dashed"),
    labels = c("predicted" = expression(
      "Analytical prediction " * frac(log(lambda), Delta*t)))
  ) +
  guides(
    colour   = guide_legend(order = 1, nrow = 1),
    linetype = guide_legend(order = 2, nrow = 1,
                            override.aes = list(colour = COL_PREDICTED))
  ) +
  labs(
    x       = expression("Observation interval  " * Delta*t),
    y       = expression("Mean bias  " * (hat(A) - A)),
    caption = sprintf(
      "%d replications per condition. Error bars: \u00B195%% CI. Lambda \u2248 %.2f throughout.",
      N_REPS_PER_COND, mean_lambda
    )
  ) +
  thesis_theme() +
  theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"))

print(p_dt)
save_thesis_plot(p_dt, "results/fig_dt.pdf",
                 width = 4.8, height = 3.2)

cat("\nFigure saved to results/fig_dt.pdf\n")