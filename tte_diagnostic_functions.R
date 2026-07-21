# Pure diagnostic helpers for binary-outcome target trial emulations.
# No data loading, database access, or study-specific paths are included.

weighted_mean <- function(x, w) sum(x * w) / sum(w)

effective_sample_size <- function(w) sum(w)^2 / sum(w^2)

overlap_weights <- function(treatment, propensity_score, epsilon = 1e-6) {
  ps <- pmin(pmax(propensity_score, epsilon), 1 - epsilon)
  ifelse(treatment == 1, 1 - ps, ps)
}

weighted_binary_effects <- function(outcome, treatment, weights) {
  ok <- !is.na(outcome) & is.finite(weights) & weights > 0
  y <- as.numeric(outcome[ok])
  z <- treatment[ok]
  w <- weights[ok]
  if (!length(y) || length(unique(z)) < 2) {
    return(data.frame(risk_treated = NA, risk_control = NA,
                      risk_difference = NA, risk_ratio = NA, odds_ratio = NA))
  }
  p1 <- weighted_mean(y[z == 1], w[z == 1])
  p0 <- weighted_mean(y[z == 0], w[z == 0])
  odds <- function(p) p / (1 - p)
  data.frame(
    risk_treated = p1,
    risk_control = p0,
    risk_difference = p1 - p0,
    risk_ratio = ifelse(p0 > 0, p1 / p0, NA_real_),
    odds_ratio = ifelse(p0 > 0 && p0 < 1 && p1 < 1, odds(p1) / odds(p0), NA_real_)
  )
}

inverse_observation_weights <- function(observation_probability,
                                        lower = 0.01, upper = 0.99) {
  probability <- pmin(pmax(observation_probability, lower), upper)
  1 / probability
}

differential_missingness_bounds <- function(outcome, treatment) {
  best_for_treated <- ifelse(is.na(outcome), 1 - treatment, outcome)
  worst_for_treated <- ifelse(is.na(outcome), treatment, outcome)
  list(best_for_treated = best_for_treated, worst_for_treated = worst_for_treated)
}

e_value <- function(risk_ratio) {
  if (!is.finite(risk_ratio) || risk_ratio <= 0) return(NA_real_)
  x <- if (risk_ratio < 1) 1 / risk_ratio else risk_ratio
  x + sqrt(x * (x - 1))
}
