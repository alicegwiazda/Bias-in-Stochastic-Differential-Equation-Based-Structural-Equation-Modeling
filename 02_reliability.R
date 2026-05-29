# =============================================================================
# 02_reliability.R
# Section 8.2: bias as a function of factor-score reliability.
#
# Two experiments that both probe the same axis of the bias formula:
#   (a) vary loadings  with sigma_eps fixed at 0.5
#   (b) vary sigma_eps with loadings  fixed at (0.8, 0.7, 0.9)
# Both reduce to changing the reliability lambda, so on a plot with lambda on
# the x-axis they trace out the same analytical curve log(lambda)/dt.
#
# Outputs:
#   results/results_reliability.rds       raw per-replication results
#   results/fig_reliability_main.pdf      headline figure: lambda on x-axis
#   results/fig_reliability_panels.pdf    supplementary: two raw knobs
# =============================================================================

source("00_helpers.R")
source("00_theme.R")

set.seed(42)

dir.create("results", showWarnings = FALSE)


# -----------------------------------------------------------------------------
# Design
# -----------------------------------------------------------------------------

N_REPS_PER_COND <- 100

LOADING_GRID <- c(0.3, 0.5, 0.7, 0.9, 1.0)   # all three indicators equal
NOISE_GRID   <- c(0.1, 0.3, 0.5, 0.7, 1.0)   # sigma_measure


# -----------------------------------------------------------------------------
# Experiment (a): vary loadings
# -----------------------------------------------------------------------------

cat("Experiment (a): varying loadings (sigma_eps = 0.5 fixed)\n")
results_loadings <- bind_rows(lapply(LOADING_GRID, function(L) {
  cat(sprintf("  loadings = %.1f\n", L))
  bind_rows(lapply(seq_len(N_REPS_PER_COND), function(i) {
    one_replication(loadings = rep(L, 3))
  }))
})) %>%
  mutate(experiment = "vary_loadings")


# -----------------------------------------------------------------------------
# Experiment (b): vary measurement noise
# -----------------------------------------------------------------------------

cat("\nExperiment (b): varying sigma_eps (loadings = (0.8,0.7,0.9) fixed)\n")
results_noise <- bind_rows(lapply(NOISE_GRID, function(s) {
  cat(sprintf("  sigma_eps = %.1f\n", s))
  bind_rows(lapply(seq_len(N_REPS_PER_COND), function(i) {
    one_replication(sigma_measure = s)
  }))
})) %>%
  mutate(experiment = "vary_noise")

results_all <- bind_rows(results_loadings, results_noise)
saveRDS(results_all, "results/results_reliability.rds")


# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("EXPERIMENT (a): bias by loading strength\n")
cat("=============================================================\n")
print(
  results_loadings %>%
    group_by(loadings_mean) %>%
    summarise(
      lambda          = round(mean(reliability),    3),
      empirical_bias  = round(mean(bias_etahat),    3),
      sd_bias         = round(sd(bias_etahat),      3),
      predicted_bias  = round(mean(predicted_bias), 3),
      n               = n()
    )
)

cat("\n=============================================================\n")
cat("EXPERIMENT (b): bias by measurement-error SD\n")
cat("=============================================================\n")
print(
  results_noise %>%
    group_by(sigma_measure) %>%
    summarise(
      lambda          = round(mean(reliability),    3),
      empirical_bias  = round(mean(bias_etahat),    3),
      sd_bias         = round(sd(bias_etahat),      3),
      predicted_bias  = round(mean(predicted_bias), 3),
      n               = n()
    )
)


# -----------------------------------------------------------------------------
# Headline figure: bias as a function of reliability (both experiments overlaid)
# -----------------------------------------------------------------------------

# Aggregate to one point per (experiment, condition). The x-axis is the
# factor-score reliability lambda, i.e. cor(eta, eta_hat)^2, which is the
# quantity the bias formula log(lambda)/dt actually refers to.
agg <- results_all %>%
  group_by(experiment, loadings_mean, sigma_measure) %>%
  summarise(
    reliability    = mean(reliability),
    mean_bias      = mean(bias_etahat),
    se_bias        = sd(bias_etahat) / sqrt(n()),
    predicted_bias = mean(predicted_bias),
    .groups        = "drop"
  )

# Build a smooth analytical curve over the empirical lambda range
lambda_range <- range(agg$reliability)
lambda_grid  <- seq(lambda_range[1], lambda_range[2], length.out = 200)
prediction_curve <- data.frame(
  reliability    = lambda_grid,
  predicted_bias = log(lambda_grid) / 1   # dt = 1
)

EXP_LABS <- c("vary_loadings" = expression("Varying loadings " * Lambda),
              "vary_noise"    = expression("Varying meas.-error SD " * sigma[epsilon]))

p_main <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted",
             linewidth = 0.4, colour = COL_REFERENCE) +
  geom_line(
    data = prediction_curve,
    aes(x = reliability, y = predicted_bias, linetype = "predicted"),
    colour = COL_PREDICTED, linewidth = 0.55
  ) +
  geom_errorbar(
    data = agg,
    aes(x = reliability, ymin = mean_bias - 1.96 * se_bias,
        ymax = mean_bias + 1.96 * se_bias),
    width = 0.012, colour = COL_ETAHAT, linewidth = 0.35
  ) +
  geom_point(
    data = agg,
    aes(x = reliability, y = mean_bias, shape = experiment),
    colour = COL_ETAHAT, size = 1.9, stroke = 0.55, fill = "white"
  ) +
  scale_shape_manual(name   = NULL,
                     values = c("vary_loadings" = 16,    # filled circle
                                "vary_noise"    = 21),   # hollow circle
                     labels = EXP_LABS) +
  scale_linetype_manual(name   = NULL,
                        values = c("predicted" = "dashed"),
                        labels = c("predicted" = expression(
                          "Analytical prediction " * frac(log(lambda), Delta*t)))) +
  guides(
    shape    = guide_legend(order = 1, nrow = 1),
    linetype = guide_legend(order = 2, nrow = 1,
                            override.aes = list(colour = COL_PREDICTED))
  ) +
  labs(
    x       = expression("Factor-score reliability  " * lambda),
    y       = expression("Mean bias  " * (hat(A) - A)),
    caption = sprintf(
      "%d replications per condition. Error bars: \u00B195%% CI. Reliability computed as cor(eta, eta-hat)^2.",
      N_REPS_PER_COND
    )
  ) +
  thesis_theme() +
  theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"))

print(p_main)
save_thesis_plot(p_main, "results/fig_reliability_main.pdf",
                 width = 4.8, height = 3.2)


# -----------------------------------------------------------------------------
# Supplementary figure: the two raw experimental knobs side by side
# -----------------------------------------------------------------------------

panel_loadings <- results_loadings %>%
  group_by(loadings_mean) %>%
  summarise(
    mean_bias_etahat = mean(bias_etahat),
    se_etahat        = sd(bias_etahat) / sqrt(n()),
    mean_bias_eta    = mean(bias_eta),
    se_eta           = sd(bias_eta) / sqrt(n()),
    .groups          = "drop"
  ) %>%
  pivot_longer(-loadings_mean,
               names_to       = c(".value", "method"),
               names_pattern  = "(mean_bias|se)_(.*)") %>%
  mutate(panel = "loadings", x = loadings_mean)

panel_noise <- results_noise %>%
  group_by(sigma_measure) %>%
  summarise(
    mean_bias_etahat = mean(bias_etahat),
    se_etahat        = sd(bias_etahat) / sqrt(n()),
    mean_bias_eta    = mean(bias_eta),
    se_eta           = sd(bias_eta) / sqrt(n()),
    .groups          = "drop"
  ) %>%
  pivot_longer(-sigma_measure,
               names_to       = c(".value", "method"),
               names_pattern  = "(mean_bias|se)_(.*)") %>%
  mutate(panel = "noise", x = sigma_measure)

panel_data <- bind_rows(panel_loadings, panel_noise) %>%
  mutate(panel = factor(panel, levels = c("loadings", "noise")))

# Labeller maps short panel codes -> plotmath expressions for strip titles
panel_labeller <- as_labeller(
  c(loadings = "Loading~strength~~Lambda",
    noise    = "Measurement-error~SD~~sigma[epsilon]"),
  default = label_parsed
)

p_panels <- ggplot(panel_data,
                   aes(x = x, y = mean_bias, colour = method, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dotted",
             linewidth = 0.4, colour = COL_REFERENCE) +
  geom_ribbon(aes(ymin = mean_bias - 1.96 * se,
                  ymax = mean_bias + 1.96 * se),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.6) +
  facet_wrap(~ panel, scales = "free_x", labeller = panel_labeller) +
  scale_colour_manual(values = c("eta" = COL_ETA, "etahat" = COL_ETAHAT),
                      labels = THESIS_LABS) +
  scale_fill_manual(values   = c("eta" = COL_ETA, "etahat" = COL_ETAHAT),
                    labels   = THESIS_LABS) +
  labs(
    x       = NULL,
    y       = expression("Mean bias  " * (hat(A) - A)),
    caption = sprintf(
      "%d replications per condition. Bands: \u00B195%% CI.",
      N_REPS_PER_COND
    )
  ) +
  thesis_theme()

print(p_panels)
save_thesis_plot(p_panels, "results/fig_reliability_panels.pdf",
                 width = 7.0, height = 3.0)

cat("\nFigures saved to results/fig_reliability_main.pdf and results/fig_reliability_panels.pdf\n")