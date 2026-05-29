# =============================================================================
# 00_helpers.R
# Shared simulation primitives used by all chapter 8 studies.
# Source this at the top of every analysis script: source("00_helpers.R")
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(psych)
})


# -----------------------------------------------------------------------------
# simulate_ou_data()
#   One replication of the data-generating process described in section 8.0.
#   Returns a tidy long-format data frame with columns:
#     id        subject identifier  (1..N)
#     time      time index          (1..Tpoints)
#     eta_true  the true latent value at this (id, time)
#     y1, y2, y3  three observed indicators
# -----------------------------------------------------------------------------
simulate_ou_data <- function(
    N             = 200,
    Tpoints       = 15,
    dt            = 1,
    A_true        = -0.5,
    sigma_process = 0.5,
    loadings      = c(0.8, 0.7, 0.9),
    sigma_measure = 0.5
) {
  
  phi <- exp(A_true * dt)        # discrete-time AR coefficient: phi = exp(A * dt)
  
  rows <- vector("list", N)
  for (i in seq_len(N)) {
    eta    <- numeric(Tpoints)
    eta[1] <- rnorm(1, 0, 1)     # initial draw from a unit-variance distribution
    for (t in 2:Tpoints) {
      eta[t] <- phi * eta[t - 1] + rnorm(1, 0, sigma_process)
    }
    rows[[i]] <- data.frame(id = i, time = seq_len(Tpoints), eta_true = eta)
  }
  data <- bind_rows(rows)
  
  # measurement model: y_k = loading_k * eta + epsilon_k
  n_rows  <- nrow(data)
  data$y1 <- loadings[1] * data$eta_true + rnorm(n_rows, 0, sigma_measure)
  data$y2 <- loadings[2] * data$eta_true + rnorm(n_rows, 0, sigma_measure)
  data$y3 <- loadings[3] * data$eta_true + rnorm(n_rows, 0, sigma_measure)
  
  data
}


# -----------------------------------------------------------------------------
# extract_factor_scores()
#   Time-naive EFA on stacked rows, with sign aligned to the true latent.
#   The sign alignment is only possible in simulation; see section 8.0.
# -----------------------------------------------------------------------------
extract_factor_scores <- function(data) {
  
  efa <- suppressWarnings(suppressMessages(
    psych::fa(
      data[, c("y1", "y2", "y3")],
      nfactors = 1,
      scores   = "regression",
      rotate   = "none"
    )
  ))
  
  eta_hat <- as.numeric(efa$scores)
  if (cor(data$eta_true, eta_hat) < 0) eta_hat <- -eta_hat
  
  data$eta_hat <- eta_hat
  data
}


# -----------------------------------------------------------------------------
# fit_ar1_drift()
#   OLS regression of y_t on y_{t-1} pooled across subjects, then
#   A_hat = log(B_hat) / dt. Returns NA if B_hat <= 0 (non-stationary draw).
# -----------------------------------------------------------------------------
fit_ar1_drift <- function(data, var, dt = 1) {
  
  d <- data %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(lagged = lag(.data[[var]])) %>%
    ungroup() %>%
    filter(!is.na(lagged))
  
  fit  <- lm(d[[var]] ~ d$lagged)
  Bhat <- coef(fit)[2]
  
  if (is.na(Bhat) || Bhat <= 0) return(NA_real_)
  log(Bhat) / dt
}


# -----------------------------------------------------------------------------
# reliability()
#   Reliability of the EFA factor scores: lambda = cor(eta, eta_hat)^2.
#   This is the lambda that appears in the bias formula log(lambda)/dt
#   derived in chapter 7. It is a property of the factor scores (the EFA
#   output), not of the indicators (the measurement-model input).
# -----------------------------------------------------------------------------
reliability <- function(data) {
  cor(data$eta_true, data$eta_hat)^2
}


# -----------------------------------------------------------------------------
# one_replication()
#   The full pipeline for one Monte Carlo draw under given conditions.
#   Returns one data.frame row; NULL if the AR fit failed.
# -----------------------------------------------------------------------------
one_replication <- function(
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
  
  lambda <- reliability(data)
  
  data.frame(
    A_true         = A_true,
    dt             = dt,
    loadings_mean  = mean(loadings),
    sigma_measure  = sigma_measure,
    reliability    = lambda,
    A_from_eta     = A_from_eta,
    A_from_etahat  = A_from_etahat,
    bias_eta       = A_from_eta    - A_true,
    bias_etahat    = A_from_etahat - A_true,
    predicted_bias = log(lambda) / dt
  )
}