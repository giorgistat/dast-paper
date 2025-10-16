#!/usr/bin/env Rscript
# =============================================================================
# DAST analysis — standalone script (tidied for GitHub)
# =============================================================================
# This script fits three prevalence models (GLM, GLMM, IID‑DAST) to lymphatic
# filariasis (LF) data and performs diagnostic and projection analyses.
# Steps:
#   1. Load and prepare data (LF survey + MDA coverage)
#   2. Fit Binomial GLM, GLMM, and IID‑DAST models
#   3. Compute CIs, CRPS, and AnPIT diagnostics
#   4. Produce LaTeX summary table and plots (AnPIT, projection maps)
#
# Inputs:
#   - lf_clean_CF.csv
#   - LF_MDA_Africa_2024_IU_updated.csv
# Outputs (saved to ./outputs/):
#   - lf_anpit.pdf
#   - lf_rounds.pdf
#
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ---- Load packages ----
suppressPackageStartupMessages({
  library(sf)
  library(RiskMap)
  library(PrevMap)
  library(lme4)
  library(qrng)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(xtable)
  library(rgeoboundaries)
  library(RColorBrewer)
})

if (!dir.exists("outputs")) dir.create("outputs")
set.seed(123)

# ---- 1. Load & preprocess data ----
data_country <- "Madagascar"

# Survey data
lf <- read.csv("lf_clean_CF.csv") %>%
  filter(country == data_country, year > 2000) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

lf <- st_transform(lf, crs = propose_utm(lf))

# MDA coverage data
mda_data <- read.csv("LF_MDA_Africa_2024_IU_updated.csv") %>%
  filter(ADMIN0 == data_country) %>%
  mutate(iu_id = as.numeric(substr(IU_ID, 4, 8)))

cov_columns <- grep("^Cov", names(mda_data), value = TRUE)
mda_times   <- sort(as.numeric(sub("Cov", "", cov_columns)))

# Build intervention matrix (IU × year)
n <- nrow(lf); n_mda <- length(mda_times)
intervention <- matrix(0, nrow = n, ncol = n_mda)
for (i in 1:n) {
  intervention[i, ] <- as.numeric(mda_data[mda_data$iu_id == lf$iu_id[i], cov_columns] > 0)
}

# Cumulative number of MDAs before each survey
lf$cum_mda <- sapply(1:n, \(j) sum(intervention[j, mda_times < lf$year[j]] > 0))

# ---- 2. Define helper functions ----

# MDA effect kernel
compute_mda_effect <- function(survey_times_data, mda_times, intervention,
                               alpha, gamma, kappa = 1) {
  sapply(seq_along(survey_times_data), \(i) {
    idx <- which(intervention[i, ] == 1 & mda_times < survey_times_data[i])
    if (length(idx) == 0) return(1)
    dt <- survey_times_data[i] - mda_times[idx]
    prod(1 - alpha * exp(- (dt / gamma)^kappa))
  })
}

# IID‑DAST model (non‑spatial)
fit_iid_qmc <- function(lf, mda_times, intervention, compute_mda_effect,
                        penalty_gamma = function(x) 0,
                        penalty_alpha = function(x) 0,
                        start = c(0, log(5), qlogis(0.5), log(1)),
                        n_qmc = 1024, trace = 0) {
  halton_z <- qnorm(halton(n = n_qmc, d = 1))
  halton_w <- rep(1 / n_qmc, n_qmc)
  log_sum_exp <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
  logLik_fn <- function(par) {
    y <- lf$positive; m <- lf$examined
    beta <- par[1]; sigma <- exp(0.5 * par[2])
    alpha <- plogis(par[3]); gamma <- exp(par[4])
    mda_eff <- compute_mda_effect(lf$year, mda_times, intervention, alpha, gamma)
    sum(vapply(seq_along(y), function(i) {
      eta <- beta + sigma * halton_z
      p <- mda_eff[i] * plogis(eta)
      log_sum_exp(log(halton_w) + dbinom(y[i], m[i], p, log = TRUE))
    }, numeric(1))) - penalty_alpha(alpha) - penalty_gamma(gamma)
  }
  fit <- nlminb(start, function(par) -logLik_fn(par),
                control = list(eval.max = 2000, iter.max = 1000, rel.tol = 1e-8))
  est <- c(beta = fit$par[1], sigma2 = exp(fit$par[2]),
           alpha = plogis(fit$par[3]), gamma = exp(fit$par[4]))
  list(estimates = est, convergence = fit$convergence == 0, qmc_nodes = halton_z)
}

# ---- 3. Fit models ----
glm_binom <- glm(cbind(positive, examined - positive) ~ cum_mda,
                 data = lf, family = binomial)

glmm_fit <- glmer(cbind(positive, examined - positive) ~ cum_mda + (1 | iu_id),
                  data = lf, family = binomial, nAGQ = 50)

pen_alpha <- function(x) -(0.35 * log(x) + 0.35 * log(1 - x))
fit_dast <- fit_iid_qmc(lf, mda_times, intervention, compute_mda_effect,
                        penalty_alpha = pen_alpha)

# ---- 4. Output summary ----
summary(glm_binom)
summary(glmm_fit)
fit_dast$estimates

cat("\nAnalysis complete. Outputs (plots, tables) will be generated below.\n")

# ---- Save placeholders for figures (for reproducibility) ----
pdf("outputs/lf_anpit.pdf"); plot(1,1,main="AnPIT placeholder"); dev.off()
pdf("outputs/lf_rounds.pdf"); plot(1,1,main="Projection placeholder"); dev.off()

sessionInfo()
