rm(list = ls())
library(sf)
library(RiskMap)
library(dplyr)
library(sp)
library(rlang)

source("simulation_functions.R")

# PARAMETERS ------------------------------------------------------------------
grid_resolution <- 0.05
points_per_region <- 25
n_sampled_per_time <- 4
time_points <- 0:3
n_times <- length(time_points)

n <- points_per_region*n_sampled_per_time*n_times

trend <- 0.5
sigma2 <- 1
phi <- 0.2

trend_means <- setNames(trend * c(-2, -1, 1, 2), c("A", "B", "C", "D"))

# SIMULATE GP & REGIONS -------------------------------------------------------
regions <- create_regions()

n_sim <- 1

res_beta <- fit_beta(0.05, 0.025, 0.95, 0.975)

penalty_f <- list(
  pn = function(x) -(res_beta$alpha*log(x) + res_beta$beta*log(1-x)),
  pn_d1 = function(x) -(res_beta$alpha/x - res_beta$beta/(1-x)),
  pn_d2 = function(x) (res_beta$alpha/x^2 + res_beta$beta/(1-x)^2)
)

grid_sim <-
  expand.grid(x = seq(0, 2, by = grid_resolution),
              y = seq(0, 2, by = grid_resolution))
n_grid <- nrow(grid_sim)

# Output
mda_grid_ext <- matrix(rep(1,5), nrow = n_grid, ncol = n_times+1)
mda_red_grid <- matrix(NA,  nrow = n_grid, ncol = n_times+1)

alpha <- 0.8
gamma <- 1/log(10*alpha/5)

mda_times <- 0:3
int_mat <- matrix(1, byrow = TRUE,
                  nrow = n_grid, ncol=n_times+1)

for(i in 1:n_grid) {
  mda_red_grid[i,] <-
    compute_mda_effect(survey_times_data = c(time_points,4),
                       mda_times = c(mda_times,4), intervention = matrix(rep(int_mat[i,],n_times+1),byrow=TRUE,
                                                                         ncol=n_times+1,
                                                                         nrow=n_times+1),
                       alpha = alpha, gamma = gamma, kappa = 1)
}

n_reg <- 4
n_par <- 8
coverage_set <- seq(0.95, 0.9999, length = 100)
n_cov <- length(coverage_set)

# Separate output lists for penalized and unpenalized models
out_pen <- out_unpen <- list()
out_pen$sampled_regions <- out_unpen$sampled_regions <- list()
out_pen$prob_true <- out_unpen$prob_true <- array(NA, dim = c(n_grid, n_times+1, n_sim))
out_pen$prob_pred <- out_unpen$prob_pred <- array(NA, dim = c(n_grid, n_times+1, n_sim))
out_pen$prob_lower <- out_unpen$prob_lower <- array(NA, dim = c(n_grid, n_times+1, n_sim, n_cov))
out_pen$prob_upper <- out_unpen$prob_upper <- array(NA, dim = c(n_grid, n_times+1, n_sim, n_cov))
out_pen$prob_true_reg <- out_unpen$prob_true_reg <- array(NA, dim = c(n_reg, n_times+1, n_sim))
out_pen$prob_pred_reg <- out_unpen$prob_pred_reg <- array(NA, dim = c(n_reg, n_times+1, n_sim))
out_pen$prob_lower_reg <- out_unpen$prob_lower_reg <- array(NA, dim = c(n_reg, n_times+1, n_sim, n_cov))
out_pen$prob_upper_reg <- out_unpen$prob_upper_reg <- array(NA, dim = c(n_reg, n_times+1, n_sim, n_cov))
out_pen$par_hat <- out_unpen$par_hat <- matrix(NA, n_par, n_sim)
out_pen$par_lower <- out_unpen$par_lower <- array(NA, dim = c(n_par, n_sim))
out_pen$par_upper <- out_unpen$par_upper <- array(NA, dim = c(n_par, n_sim))

reg_names <- c("A", "B", "C", "D")
list_reg <- sapply(reg_names, function(x) which(assign_regions_to_grid(grid_sim)[,3]==x))
time_pred_set <- 0:4

for(i in 1:n_sim) {

  gp_obj <- simulate_gp(
    grid_res = grid_resolution,
    sigma2 = sigma2,
    phi = phi,
    regions = regions,
    trend_means = trend_means
  )

  gp_grid_base <- gp_obj$gp_grid_base
  gp_grid_with_trend <- gp_obj$gp_grid_with_trend
  coords_mat <- gp_obj$coords_mat
  unit_grid <- gp_obj$unit_grid

  final_samples <- simulate_sampling(
    regions = regions,
    time_points = time_points,
    n_sampled = n_sampled_per_time,
    grid_res = grid_resolution,
    n_grid_points = points_per_region,
    unit_grid_df = unit_grid,
    coords_mat = coords_mat,
    trend_means = trend_means)

  sampled_regions <- data.frame(matrix(0,length(trend_means), n_times))
  colnames(sampled_regions) <- 0:3
  rownames(sampled_regions) <- c("A","B","C","D")

  for(j in 1:n_times) {
    reg_j <- as.character(unique(final_samples[final_samples$time==time_points[j],]$region))
    sampled_regions[reg_j,j] <- 1
  }

  for(model_type in c("unpen", "pen")) {
    use_penalty <- (model_type == "pen")
    samples <- final_samples
    samples$mda_red <- compute_mda_effect(survey_times_data = samples$time,
                                          mda_times = mda_times,
                                          intervention = int_mat,
                                          alpha = alpha, gamma = gamma, kappa = 1)
    samples$prob_star_value <- 1 / (1 + exp(-samples$lp_value))
    samples$prob_value <- samples$prob_star_value * samples$mda_red
    samples$den <- 50
    samples$resp <- rbinom(n, size = samples$den, prob = samples$prob_value)
    samples$region <- droplevels(factor(samples$region, levels = c("A", "B", "C", "D")))

    par0 <- list(beta = trend_means[as.character(sort(unique(samples$region)))],
                 sigma2 = sigma2, phi = phi, alpha = alpha, gamma = gamma)

    dast_fit <- dast(resp ~ -1 + region + gp(x,y,time), crs = 3857, scale_to_km = FALSE,
                     den = den, par0 = par0, start_pars = par0,
                     mda_times = mda_times,
                     int_mat = int_mat, power_val = 1,
                     data = st_drop_geometry(samples), messages = FALSE,
                     penalty = if(use_penalty) penalty_f else NULL)

    out_list <- if(use_penalty) out_pen else out_unpen

    levels(samples$region)

    sel_reg_i <- sapply(reg_names, function(x) any(x==levels(samples$region)))
    sel_par_i <- c(sel_reg_i,rep(TRUE, length(coef(dast_fit))-1))
    out_list$par_hat[sel_par_i,i] <- unlist(coef(dast_fit))
    sum_dast <- summary(dast_fit)
    out_list$par_lower[sel_par_i,i] <- c(sum_dast$reg_coef[,2], sum_dast$sp[,2], sum_dast$dast_par[,2])
    out_list$par_upper[sel_par_i,i] <- c(sum_dast$reg_coef[,3], sum_dast$sp[,3], sum_dast$dast_par[,3])

    prob_star_grid <- 1 / (1 + exp(-gp_obj$gp_grid_with_trend$gp_with_trend))
    out_list$prob_true[,,i] <- as.matrix(mda_red_grid) * prob_star_grid
    out_list$prob_true_reg[,,i] <- t(sapply(1:n_reg, function(j)
      apply(out_list$prob_true[list_reg[[j]],,i], 2, mean)))

    grid_pred <- st_coordinates(gp_grid_base)
    colnames(grid_pred) <- c("x","y")
    counts_reg <-   apply(sampled_regions,1,sum)
    no_reg <- any(counts_reg==0)
    predictors_grid <- data.frame(region = factor(gp_grid_with_trend$region,
                                                  levels = reg_names))

    if(no_reg) {
      which_no_reg <- names(which(counts_reg==0))
      inc <- sapply(gp_grid_with_trend$region, function(x) prod(x!=which_no_reg)==1)
      grid_pred <- grid_pred[inc,]
      predictors_grid <- data.frame(region = factor(predictors_grid[inc,],
                                                    levels = c("A", "B", "C", "D")))
      predictors_grid$region <- droplevels(predictors_grid$region)

    } else {
      inc <- rep(TRUE, n_grid)
    }
    inc_reg <- which(apply(sampled_regions,1,sum)>0)


    pred_dast_m <-
      pred_over_grid(dast_fit, grid_pred = st_as_sfc(gp_grid_base[inc,]),
                     predictors = predictors_grid,
                     type = "joint", messages = FALSE)

    for(h in 1:(n_times+1)) {
      pred_res <-
        pred_target_grid(pred_dast_m, time_pred = time_pred_set[h],
                         f_target = list(prev = function(x) 1/(1+exp(-x))),
                         pd_summary = list(mean=mean,
                                           lower = function(x) quantile(x, (1-coverage_set)/2),
                                           upper = function(x) quantile(x, 1-(1-coverage_set)/2)),
                         mda_grid = mda_grid_ext[inc,])
      out_list$prob_pred[inc,h,i] <- pred_res$target$prev$mean
      out_list$prob_lower[inc,h,i,] <- pred_res$target$prev$lower
      out_list$prob_upper[inc,h,i,] <- pred_res$target$prev$upper

      pred_reg_res <- pred_target_shp(pred_dast_m,
                                      shp = district_sf[inc_reg,],
                                      time_pred = time_pred_set[h],
                                      f_target = list(prev = function(x) 1/(1+exp(-x))),
                                      standardize_weights = FALSE,
                                      pd_summary = list(mean=mean,
                                                        lower = function(x) quantile(x, (1-coverage_set)/2),
                                                        upper = function(x) quantile(x, 1-(1-coverage_set)/2)),
                                      mda_grid = mda_grid_ext[inc,],
                                      return_shp = FALSE,
                                      messages = FALSE)
      n_inc_reg <- length(inc_reg)
      out_list$prob_pred_reg[inc_reg,h,i] <- sapply(1:n_inc_reg, function(k) pred_reg_res$target[[k]]$prev$mean)
      out_list$prob_lower_reg[inc_reg,h,i,] <- sapply(1:n_inc_reg, function(k) pred_reg_res$target[[k]]$prev$lower)
      out_list$prob_upper_reg[inc_reg,h,i,] <- sapply(1:n_inc_reg, function(k) pred_reg_res$target[[k]]$prev$upper)
    }

    if(use_penalty) {
      out_pen <- out_list
    } else {
      out_unpen <- out_list
    }
  }

  cat("Iter:", i, "completed\n")
}

save(out_pen, out_unpen, file = paste("sim1_", floor(10e5*runif(1)), ".RData", sep=""))
