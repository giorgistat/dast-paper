###############################################################################
# LF DAST / GLMM / Binomial analysis + projections and MDA-rounds map
#
# - Loads LF data and MDA coverage data for Madagascar
# - Constructs intervention matrix & cumulative MDA covariate
# - Fits:
#     (i)   Binomial GLM
#     (ii)  Binomial GLMM with random intercepts (id_loc)
#     (iii) IID–DAST model (QMC integration)
# - Computes:
#     * Parametric bootstrap CIs for DAST parameters
#     * Profile-likelihood CI for GLMM random-effect variance σ²
#     * Cross-validated CRPS for all three models
#     * Average non-randomised PIT (AnPIT) curves (saved to PDF)
#     * LaTeX comparison table (printed to console)
# - Uses a previously fitted DAST model on a grid (dast_m) to:
#     * Project prevalence under alternative MDA schedules
#     * Derive number of additional MDA rounds needed to reach <1% prevalence
#     * Plot the resulting “rounds needed” map for Madagascar (saved to PDF)
#
# Dependencies: RiskMap, PrevMap, lme4, qrng, sf, dplyr, tidyr, purrr, ggplot2,
#              xtable, rgeoboundaries, jsonlite
# External files (expected in working directory):
#   - "lf_clean_CF.csv"
#   - "LF_MDA_Africa_2024_IU_updated.csv"
#   - "simulation_functions.R"     (must define pred_over_grid(), pred_target_grid(), etc.)
#   - "mgd_fit.RData"              (must contain object `dast_m`)
###############################################################################

# ----------------------------- 0. Workspace ---------------------------------
rm(list = ls())

# CRAN / GitHub packages ------------------------------------------------------
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
library(jsonlite)

set.seed(123)

# --------------------------- 1. Load & wrangle data --------------------------

## choose country -------------------------------------------------------------
data_country <- "Madagascar"

## survey prevalence data -----------------------------------------------------
lf <- read.csv("lf_clean_CF.csv") %>%
  filter(country == data_country) %>%
  filter(year > 2000) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

## project to appropriate UTM -------------------------------------------------
data_utm <- propose_utm(lf)
lf <- st_transform(lf, crs = data_utm)

## MDA administrative coverage data ------------------------------------------
mda_data <- read.csv("LF_MDA_Africa_2024_IU_updated.csv") %>%
  filter(ADMIN0 == data_country) %>%
  mutate(iu_id = as.numeric(substr(IU_ID, 4, 8)))

## identify coverage columns & years -----------------------------------------
mda_columns <- grepl("^Cov", names(mda_data))
cov_columns <- names(mda_data)[mda_columns]
mda_times   <- sort(as.numeric(sub("Cov", "", cov_columns)))

## intervention matrix (IU × MDA year) ---------------------------------------
n_mda <- length(mda_times)
n     <- nrow(lf)
intervention <- matrix(NA, ncol = n_mda, nrow = n)
for (i in 1:n) {
  intervention[i, ] <-
    1 * (as.numeric(mda_data[which(mda_data$iu_id == lf$iu_id[i]), mda_columns]) > 0)
}

## cumulative number of MDAs before each survey ------------------------------
lf$cum_mda <- sapply(
  1:nrow(lf),
  function(j) sum(intervention[j, mda_times < lf$year[j]] > 0)
)

## empirical logit ------------------------------------------------------------
lf$elogit <- with(lf, log((positive + 0.5) / (examined - positive + 0.5)))

## coordinates + random-effects ID for GLMM ----------------------------------
lf <- lf %>%
  mutate(
    x = st_coordinates(lf)[, 1],
    y = st_coordinates(lf)[, 2]
  )
lf$id_loc <- create.ID.coords(st_drop_geometry(lf), ~ x + y)

# ---------------------- 2. Helper functions (all inline) --------------------

## 2.1  effect of prior MDAs --------------------------------------------------
compute_mda_effect <- function(survey_times_data, mda_times, intervention,
                               alpha, gamma, kappa = 1) {
  out <- rep(1, length(survey_times_data))
  for (i in seq_along(survey_times_data)) {
    idx <- which(intervention[i, ] == 1 & mda_times < survey_times_data[i])
    if (length(idx) > 0) {
      dt <- survey_times_data[i] - mda_times[idx]
      out[i] <- prod(1 - alpha * exp(- (dt / gamma)^kappa))
    }
  }
  out
}

## 2.2  IID-DAST log-likelihood with Sobol/Halton QMC ------------------------
fit_iid_qmc <- function(lf, mda_times, intervention,
                        compute_mda_effect,
                        penalty_gamma = function(x) 0,
                        penalty_alpha = function(x) 0,
                        start = c(0, log(5), qlogis(0.5), log(1)),
                        n_qmc = 1024,
                        trace = 1) {
  
  ## --- Extract location IDs ---
  id_loc <- lf$id_loc
  unique_locs <- unique(id_loc)
  n_locs <- length(unique_locs)
  
  ## --- Halton QMC nodes (1D) ---
  halton_z <- qnorm(halton(n = n_qmc, d = 1))
  halton_w <- rep(1 / n_qmc, n_qmc)
  
  ## --- Stable log-sum-exp ---
  log_sum_exp <- function(v) {
    m <- max(v)
    m + log(sum(exp(v - m)))
  }
  
  ## --- Log-likelihood function ---
  logLik_fn <- function(par) {
    y        <- lf$positive
    units_m  <- lf$examined
    n        <- length(y)
    
    beta     <- par[1]
    sigma    <- exp(0.5 * par[2])
    alpha    <- plogis(par[3])
    gamma    <- exp(par[4])          # decay time
    
    mu_fix   <- rep(beta, n)
    mda_eff  <- compute_mda_effect(
      survey_times_data = lf$year,
      mda_times         = mda_times,
      intervention      = intervention,
      alpha             = alpha,
      gamma             = gamma,
      kappa             = 1
    )
    
    ## --- Compute log-likelihood by location ---
    logLik <- sum(
      vapply(seq_along(unique_locs), function(loc_idx) {
        # Get indices for this location
        loc_id <- unique_locs[loc_idx]
        idx_loc <- which(id_loc == loc_id)
        
        # For each QMC node (random effect value)
        log_fk <- vapply(seq_along(halton_z), function(k) {
          u_k <- halton_z[k]
          
          # Contribution from all observations at this location
          log_lik_loc <- sum(
            vapply(idx_loc, function(i) {
              eta_i <- mu_fix[i] + sigma * u_k
              p_i   <- mda_eff[i] * plogis(eta_i)
              dbinom(y[i], size = units_m[i], prob = p_i, log = TRUE)
            }, numeric(1))
          )
          
          log(halton_w[k]) + log_lik_loc
        }, numeric(1))
        
        # Integrate over random effect for this location
        log_sum_exp(log_fk)
      }, numeric(1))
    )
    
    ## --- total penalised log-likelihood ---
    logLik - penalty_alpha(alpha) - penalty_gamma(gamma)
  }
  
  ## --- Run optimisation ---
  raw_fit <- nlminb(
    start,
    objective = function(par) -logLik_fn(par),
    control   = list(trace = trace, eval.max = 2000,
                     iter.max = 1000, rel.tol = 1e-8)
  )
  
  ## --- Transform parameters to natural scale ---
  par_opt <- raw_fit$par
  est <- c(
    beta    = par_opt[1],
    sigma2  = exp(par_opt[2]),       # variance scale
    alpha   = plogis(par_opt[3]),
    gamma   = exp(par_opt[4])
  )
  
  ## --- Return structured output ---
  list(
    estimates   = est,
    objective   = raw_fit$objective,
    converged   = raw_fit$convergence == 0,
    raw_optim   = raw_fit,
    logLik_fn   = logLik_fn,
    qmc_nodes   = halton_z
  )
}


## 2.3  CRPS utility ----------------------------------------------------------
crps_sample <- function(y, x) {
  mean(abs(x - y)) - 0.5 * mean(abs(outer(x, x, "-")))
}

## 2.4  Yearly cross-validated CRPS for 3 models -----------------------------
compute_crps_all <- function(lf, intervention, mda_times,
                             glm_binom, glmm_fit, dast_fit,
                             compute_mda_effect, prop_test = 0.2, M = 400) {
  set.seed(99)
  lf$test_set <- FALSE
  for (yr in unique(lf$year)) {
    idx <- which(lf$year == yr)
    lf$test_set[sample(idx, ceiling(prop_test * length(idx)))] <- TRUE
  }
  idx    <- which(lf$test_set)
  y_true <- lf$positive[idx]
  n_obs  <- lf$examined[idx]
  
  # Binomial draws -----------------------------------------------------------
  p_bin <- plogis(predict(glm_binom, lf[idx, ], type = "link"))
  y_bin <- replicate(M, rbinom(length(idx), n_obs, p_bin))
  
  # GLMM draws ---------------------------------------------------------------
  b0  <- fixef(glmm_fit)[1]
  b1  <- fixef(glmm_fit)["cum_mda"]
  sdU <- sqrt(as.data.frame(VarCorr(glmm_fit))$vcov[1])
  u   <- matrix(rnorm(M * length(idx), 0, sdU), nrow = M)
  p_glm <- plogis(sweep(u, 2, b0 + b1 * lf$cum_mda[idx], "+"))
  y_glm <- matrix(
    rbinom(M * length(idx), rep(n_obs, each = M), c(p_glm)),
    nrow = M
  )
  
  # DAST draws ---------------------------------------------------------------
  beta  <- dast_fit$estimates["beta"]
  sigma <- sqrt(dast_fit$estimates["sigma2"])
  alpha <- dast_fit$estimates["alpha"]
  gamma <- dast_fit$estimates["gamma"]
  z     <- dast_fit$qmc_nodes
  w     <- rep(1 / length(z), length(z))
  
  mda_eff <- compute_mda_effect(
    lf$year[idx], mda_times, intervention[idx, ], alpha, gamma
  )
  
  z_samp <- matrix(sample(z, M * length(idx), TRUE, prob = w), nrow = M)
  p_da   <- sweep(plogis(beta + sigma * z_samp), 2, mda_eff, `*`)
  y_da   <- matrix(
    rbinom(M * length(idx), rep(n_obs, each = M), c(p_da)),
    nrow = M
  )
  
  CRPS <- tibble(
    Binomial = vapply(seq_along(idx), function(j)
      crps_sample(y_true[j], y_bin[, j]), numeric(1)),
    GLMM     = vapply(seq_along(idx), function(j)
      crps_sample(y_true[j], y_glm[, j]), numeric(1)),
    DAST     = vapply(seq_along(idx), function(j)
      crps_sample(y_true[j], y_da[, j]), numeric(1))
  )
  
  colMeans(CRPS)
}

# --------------------------- 3. Fit three models -----------------------------

## 3.1  Standard binomial GLM -------------------------------------------------
glm_binom <- glm(
  cbind(positive, examined - positive) ~ cum_mda,
  data = lf, family = binomial
)

## 3.2  GLMM ------------------------------------------------------------------
glmm_fit <- glmer(
  cbind(positive, examined - positive) ~ cum_mda + (1 | id_loc),
  data = lf, family = binomial, nAGQ = 100
)

# σ² estimate & profile-likelihood CI ----------------------------------------
ci_all <- suppressMessages(confint(glmm_fit, method = "profile"))
sd_ci  <- ci_all[grep(".sig01", rownames(ci_all)), , drop = FALSE]
varU_hat <- as.data.frame(VarCorr(glmm_fit))$vcov[1]
varU_ci  <- sd_ci^2

## 3.3  IID-DAST --------------------------------------------------------------
pen_alpha <- function(x) -(0.35 * log(x) + 0.35 * log(1 - x))

fit_dast <- fit_iid_qmc(
  lf, mda_times, intervention, compute_mda_effect,
  penalty_alpha = pen_alpha,
  penalty_gamma = function(x) 0,
  trace = 1
)

# --------------------------- 4. Parametric bootstrap ------------------------
set.seed(202)
B <- 1000

boot_est <- replicate(B, {
  beta  <- fit_dast$estimates["beta"]
  sigma <- sqrt(fit_dast$estimates["sigma2"])
  alpha <- fit_dast$estimates["alpha"]
  gamma <- fit_dast$estimates["gamma"]
  
  z <- sigma * rnorm(nrow(lf))
  mda_eff <- compute_mda_effect(lf$year, mda_times, intervention, alpha, gamma)
  p_draws <- plogis(beta + z) * mda_eff
  
  y_new <- rbinom(nrow(lf), lf$examined, p_draws)
  lf_boot <- lf
  lf_boot$positive <- y_new
  
  fit_b <- tryCatch(
    fit_iid_qmc(
      lf_boot, mda_times, intervention, compute_mda_effect,
      penalty_alpha = pen_alpha,
      penalty_gamma = function(x) 0,
      trace = 0
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit_b$converged) && fit_b$converged) {
    return(unlist(fit_b$estimates))
  } else {
    return(rep(NA, 4))
  }
}, simplify = "matrix")

boot_df  <- as.data.frame(t(boot_est))
boot_cis <- t(apply(boot_df, 2, quantile,
                    probs = c(0.025, 0.975), na.rm = TRUE))

# --------------------------- 5. Metrics summary -----------------------------
crps_all <- compute_crps_all(
  lf, intervention, mda_times,
  glm_binom, glmm_fit, fit_dast,
  compute_mda_effect
)

# ---- 5.1 Extract estimates and CIs -----------------------------------------

## GLM (Binomial)
beta_binom <- coef(glm_binom)["(Intercept)"]
cum_binom  <- coef(glm_binom)["cum_mda"]
ci_glm     <- confint(glm_binom)

## GLMM
beta_glmm   <- fixef(glmm_fit)["(Intercept)"]
cum_glmm    <- fixef(glmm_fit)["cum_mda"]
ci_glmm_fix <- confint(glmm_fit,
                       parm = c("(Intercept)", "cum_mda"),
                       method = "Wald")
ci_beta_glmm <- ci_glmm_fix["(Intercept)", ]
ci_cum_glmm  <- ci_glmm_fix["cum_mda", ]

varU_glmm <- as.numeric(VarCorr(glmm_fit)[["id_loc"]])
ci_varU   <- varU_ci  # from profile on SD, squared to variance

## DAST
beta_dast   <- fit_dast$estimates["beta"]
alpha_dast  <- fit_dast$estimates["alpha"]
gamma_dast  <- fit_dast$estimates["gamma"]
sigma2_dast <- fit_dast$estimates["sigma2"]
ci_nat      <- boot_cis

# ---- 5.2 Assemble LaTeX-ready comparison table -----------------------------

latex_table <- data.frame(
  Parameter = c("beta", "cum\\_mda", "alpha", "gamma", "Variance", "CRPS"),
  
  Binomial = c(
    sprintf("%.3f (%.3f–%.3f)",
            beta_binom, ci_glm["(Intercept)", 1], ci_glm["(Intercept)", 2]),
    sprintf("%.3f (%.3f–%.3f)",
            cum_binom,  ci_glm["cum_mda", 1], ci_glm["cum_mda", 2]),
    "--",
    "--",
    "--",
    sprintf("%.3f", crps_all["Binomial"])
  ),
  
  GLMM = c(
    sprintf("%.3f (%.3f–%.3f)",
            beta_glmm, ci_glmm_fix[1, 1], ci_glmm_fix[1, 2]),
    sprintf("%.3f (%.3f–%.3f)",
            cum_glmm, ci_glmm_fix[2, 1], ci_glmm_fix[2, 2]),
    "--",
    "--",
    sprintf("%.3f (%.3f–%.3f)",
            varU_glmm, ci_varU[1], ci_varU[2]),
    sprintf("%.3f", crps_all["GLMM"])
  ),
  
  DAST = c(
    sprintf("%.3f (%.3f–%.3f)",
            beta_dast, ci_nat["beta", 1],  ci_nat["beta", 2]),
    "--",
    sprintf("%.3f (%.3f–%.3f)",
            alpha_dast, ci_nat["alpha", 1], ci_nat["alpha", 2]),
    sprintf("%.3f (%.3f–%.3f)",
            gamma_dast, ci_nat["gamma", 1], ci_nat["gamma", 2]),
    sprintf("%.3f (%.3f–%.3f)",
            sigma2_dast, ci_nat["sigma2", 1], ci_nat["sigma2", 2]),
    sprintf("%.3f", crps_all["DAST"])
  ),
  stringsAsFactors = FALSE
)

print(
  xtable(latex_table, align = "llccc", digits = 3),
  include.rownames = FALSE,
  sanitize.text.function = identity
)

# ----------------------- 6. AnPIT computation & plot ------------------------

# 6.1  Define hold-out indices (same rule as compute_crps_all) ---------------
prop_test <- 0.10
lf$test_set <- FALSE
for (yr in unique(lf$year)) {
  idx_yr <- which(lf$year == yr)
  lf$test_set[sample(idx_yr, ceiling(prop_test * length(idx_yr)))] <- TRUE
}
test_idx <- lf$test_set

u_grid <- seq(0, 1, length.out = 1000)

# 6.2  Helper: GLM AnPIT -----------------------------------------------------
compute_anpit_glm <- function(glm_fit, lf, test_idx,
                              u_grid = seq(0, 1, length.out = 1000)) {
  idx      <- which(test_idx)
  y_obs    <- lf$positive[idx]
  n_obs    <- lf$examined[idx]
  p_hat    <- plogis(predict(glm_fit, newdata = lf[idx, ], type = "link"))
  
  AnPIT_mat <- sapply(seq_along(idx), function(i) {
    Fy_minus <- pbinom(y_obs[i] - 1L, size = n_obs[i], prob = p_hat[i])
    Fy       <- pbinom(y_obs[i],       size = n_obs[i], prob = p_hat[i])
    ifelse(
      u_grid <= Fy_minus, 0,
      ifelse(
        u_grid <= Fy,
        (u_grid - Fy_minus) /
          pmax(Fy - Fy_minus, .Machine$double.eps),
        1
      )
    )
  })
  
  data.frame(
    u     = u_grid,
    AnPIT = rowMeans(AnPIT_mat, na.rm = TRUE),
    Model = "GLM"
  )
}

# 6.3  Helper: GLMM AnPIT ----------------------------------------------------
compute_anpit_glmm <- function(glmer_fit, lf, test_idx,
                               u_grid = seq(0, 1, length.out = 1000),
                               n_qmc = 300) {
  b       <- fixef(glmer_fit)
  beta0   <- b[1]
  beta1   <- if ("cum_mda" %in% names(b)) b["cum_mda"] else b[2]
  sigma_u <- sqrt(as.data.frame(VarCorr(glmer_fit))$vcov[1])
  z <- qnorm(halton(n = n_qmc, d = 1))
  w <- rep(1 / n_qmc, n_qmc)
  
  idx      <- which(test_idx)
  y_obs    <- lf$positive[idx]
  n_obs    <- lf$examined[idx]
  cum_mda  <- lf$cum_mda[idx]
  
  mixture_cdf <- function(y, i) {
    linpred <- beta0 + beta1 * cum_mda[i] + sigma_u * z
    p_k     <- plogis(linpred)
    sum(w * pbinom(y, size = n_obs[i], prob = p_k))
  }
  
  AnPIT_mat <- sapply(seq_along(idx), function(i) {
    Fy_minus <- mixture_cdf(y_obs[i] - 1L, i)
    Fy       <- mixture_cdf(y_obs[i],       i)
    ifelse(
      u_grid <= Fy_minus, 0,
      ifelse(
        u_grid <= Fy,
        (u_grid - Fy_minus) /
          pmax(Fy - Fy_minus, .Machine$double.eps),
        1
      )
    )
  })
  
  data.frame(
    u     = u_grid,
    AnPIT = rowMeans(AnPIT_mat, na.rm = TRUE),
    Model = "GLMM"
  )
}

# 6.4  Helper: DAST AnPIT ----------------------------------------------------
compute_anpit_dast <- function(dast_fit, lf, test_idx, int_mat,
                               mda_times, compute_mda_effect,
                               u_grid = seq(0, 1, length.out = 1000),
                               n_qmc = 300) {
  beta  <- dast_fit$estimates["beta"]
  sigma <- sqrt(dast_fit$estimates["sigma2"])
  alpha <- dast_fit$estimates["alpha"]
  gamma <- dast_fit$estimates["gamma"]
  
  z <- qnorm(halton(n = n_qmc, d = 1))
  w <- rep(1 / n_qmc, n_qmc)
  
  idx      <- which(test_idx)
  y_obs    <- lf$positive[idx]
  n_obs    <- lf$examined[idx]
  years    <- lf$year[idx]
  
  mda_eff <- compute_mda_effect(
    years, mda_times,
    int_mat[idx, ], alpha, gamma, kappa = 1
  )
  
  mixture_cdf <- function(y, i) {
    p_k <- mda_eff[i] * plogis(beta + sigma * z)
    sum(w * pbinom(y, size = n_obs[i], prob = p_k))
  }
  
  AnPIT_mat <- sapply(seq_along(idx), function(i) {
    Fy_minus <- mixture_cdf(y_obs[i] - 1L, i)
    Fy       <- mixture_cdf(y_obs[i],       i)
    ifelse(
      u_grid <= Fy_minus, 0,
      ifelse(
        u_grid <= Fy,
        (u_grid - Fy_minus) /
          pmax(Fy - Fy_minus, .Machine$double.eps),
        1
      )
    )
  })
  
  data.frame(
    u     = u_grid,
    AnPIT = rowMeans(AnPIT_mat, na.rm = TRUE),
    Model = "DAST"
  )
}

# 6.5  Compute AnPIT curves ---------------------------------------------------
anpit_glm  <- compute_anpit_glm(glm_binom, lf, test_idx)
anpit_glmm <- compute_anpit_glmm(glmm_fit, lf, test_idx)
anpit_dast <- compute_anpit_dast(
  fit_dast, lf, test_idx,
  intervention, mda_times, compute_mda_effect
)

anpit_plot_df <- bind_rows(anpit_glm, anpit_glmm, anpit_dast)

# 6.6  Plot AnPIT to PDF ------------------------------------------------------
pdf("lf_anpit.pdf")
anpit_plot <- ggplot(anpit_plot_df, aes(u, AnPIT, colour = Model)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0,
              colour = "black", linetype = "dashed") +
  coord_equal() +
  labs(
    x = expression(u),
    y = "Average non-randomised Probability Integral Transform (AnPIT)"
  ) +
  scale_colour_manual(values = c(
    "GLM"  = "forestgreen",
    "GLMM" = "blue",
    "DAST" = "red"
  )) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
print(anpit_plot)
dev.off()

# ---------------- 7. Projections & MDA rounds needed (<1% prev)  ------------

# For projections we use a previously fitted DAST model on a prediction grid
# stored in "mgd_fit.RData" (object `dast_m`) and helpers in "simulation_functions.R".

source("simulation_functions.R")
load("mgd_fit.RData")   # should create object `dast_m`

pred_dast_loc <- pred_over_grid(dast_m)

year_set <- 2023:2025

# Create list of alternative future MDA patterns (0,1,2,3 extra annual rounds)
intervention_list <- lapply(0:3, function(i) {
  cbind(
    intervention,
    matrix(
      c(rep(1, i), rep(0, 4 - i)),
      byrow = TRUE,
      ncol  = 4,
      nrow  = nrow(intervention)
    )
  )
})

n_int   <- 4
n_year  <- length(year_set)
n_loc   <- nrow(lf)
pred_out <- array(NA, dim = c(n_loc, n_year, n_int))

# extend mda_times to include projection years -------------------------------
pred_dast_loc$mda_times <- c(pred_dast_loc$mda_times, year_set)

for (i in 1:n_year) {
  for (j in 1:n_int) {
    out_ij <- pred_target_grid(
      pred_dast_loc,
      time_pred = year_set[i],
      mda_grid  = intervention_list[[j]],
      f_target  = list(prev = function(x) 1 / (1 + exp(-x)))
    )
    pred_out[, i, j] <- out_ij$target$prev$mean
  }
}

# pred_out: 3D array of predicted prevalence [location, year, intervention]
# x: threshold reduction target (e.g., 0.01 = 1% prevalence)
# baseline_idx: index of year to use as baseline (e.g., 1 for first year)

check_threshold_reduction <- function(pred_out, x, baseline_idx = 1) {
  n <- dim(pred_out)[1]  # number of locations
  T <- dim(pred_out)[2]  # number of time points
  J <- dim(pred_out)[3]  # number of interventions
  
  results <- character(n)
  
  for (i in 1:n) {
    baseline <- pred_out[i, baseline_idx, 1]
    if (baseline < x) {
      results[i] <- "Already below"
      next
    }
    
    found <- FALSE
    for (j in 1:J) {
      if (any(pred_out[i, 1:3, j] < x)) {
        results[i] <- paste0(j, " rounds")
        found <- TRUE
        break
      }
    }
    
    if (!found) {
      results[i] <- "> 3 rounds"
    }
  }
  
  results
}

lf$mda_rounds_predicted <- check_threshold_reduction(pred_out, x = 0.01)

# ---------------- 8. Map of predicted MDA rounds (Madagascar) ---------------

# ADM0 boundary from GeoBoundaries (MDG)
api_url <- "https://www.geoboundaries.org/api/current/gbOpen/MDG/ADM0/"
meta         <- jsonlite::fromJSON(api_url)
mada_admin0  <- st_read(meta$gjDownloadURL, quiet = TRUE)
mada_admin0  <- st_transform(mada_admin0, crs = st_crs(lf))

# Ensure factor includes all labels we want to display -----------------------
lf$mda_rounds_predicted <- factor(
  lf$mda_rounds_predicted,
  levels = c("Already below", "1 rounds", "2 rounds", "3 rounds", "> 3 rounds")
)

# Optionally relabel "1 rounds" → "1 round" for nicer legend text
levels(lf$mda_rounds_predicted)[levels(lf$mda_rounds_predicted) == "1 rounds"] <- "1 round"

lf$mda_rounds_predicted <- factor(
  lf$mda_rounds_predicted,
  levels = c("Already below", "1 round", "2 rounds", "3 rounds", "> 3 rounds")
)

# One unique colour per label -----------------------------------------------
fill_colors <- c(
  "Already below" = "grey70",
  "1 round"       = "lightgreen",
  "2 rounds"      = "skyblue",
  "3 rounds"      = "dodgerblue",
  "> 3 rounds"    = "red"
)

# Extract coordinates for plotting with geom_point ---------------------------
lf_points <- lf %>%
  mutate(
    x = st_coordinates(lf)[, 1],
    y = st_coordinates(lf)[, 2]
  ) %>%
  st_drop_geometry()

# Plot and save to PDF -------------------------------------------------------
pdf("predicted_rounds.pdf")
ggplot() +
  geom_sf(data = mada_admin0, fill = NA,
          colour = "black", linewidth = 0.6) +
  
  geom_point(
    data   = lf_points,
    aes(x = x, y = y, fill = mda_rounds_predicted),
    shape  = 21, size = 2,
    colour = "white", stroke = 0.2
  ) +
  
  # Invisible point to ensure "1 round" appears in legend (robust to empty class)
  geom_point(
    data = data.frame(
      x = Inf, y = Inf,
      mda_rounds_predicted = factor(
        "1 round",
        levels = c("Already below", "1 round", "2 rounds", "3 rounds", "> 3 rounds")
      )
    ),
    aes(x, y, fill = mda_rounds_predicted),
    shape = 21, size = 0,
    show.legend = TRUE,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = fill_colors,
    name   = "Rounds needed",
    drop   = FALSE,
    guide  = guide_legend(
      override.aes = list(shape = 21, colour = "white", size = 4)
    )
  ) +
  coord_sf(crs = st_crs(lf)) +
  theme_minimal() +
  labs(
    title    = "Predicted MDA rounds needed to reduce LF prevalence <1%",
    subtitle = "5-year projection under alternative MDA schedules",
    caption  = "ADM0 boundary from GeoBoundaries"
  )
dev.off()

# quick summary of mean baseline prevalence by predicted category ------------
lf$prev <- lf$positive / lf$examined
print(tapply(lf$prev, lf$mda_rounds_predicted, mean, na.rm = TRUE))
