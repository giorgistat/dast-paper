# simulation_functions.R

library(sf)
library(dplyr)
library(MASS)
library(RANN)

# Your coordinates
district_coords <- list(
  A = matrix(c(0, 1, 1, 0, 0, 0, 0, 1, 1, 0), ncol = 2),
  B = matrix(c(1, 2, 2, 1, 1, 0, 0, 1, 1, 0), ncol = 2),
  C = matrix(c(0, 1, 1, 0, 0, 1, 1, 2, 2, 1), ncol = 2),
  D = matrix(c(1, 2, 2, 1, 1, 1, 1, 2, 2, 1), ncol = 2)
)

# Convert to list of POLYGON geometries
polygons <- lapply(names(district_coords), function(name) {
  coords <- district_coords[[name]]
  st_polygon(list(coords))
})

# Create sf object
district_sf <- st_sf(
  district = names(district_coords),
  geometry = st_sfc(polygons)
)

# Set CRS if needed, e.g., WGS84
st_crs(district_sf) <- 3857

# Create square regions A–D
create_regions <- function() {
  district_coords <- list(
    A = matrix(c(0, 1, 1, 0, 0, 0, 0, 1, 1, 0), ncol = 2),
    B = matrix(c(1, 2, 2, 1, 1, 0, 0, 1, 1, 0), ncol = 2),
    C = matrix(c(0, 1, 1, 0, 0, 1, 1, 2, 2, 1), ncol = 2),
    D = matrix(c(1, 2, 2, 1, 1, 1, 1, 2, 2, 1), ncol = 2)
  )

  regions_sf <- lapply(names(district_coords), function(name) {
    poly <- st_polygon(list(district_coords[[name]]))
    st_sf(name = name, geometry = st_sfc(poly), crs = 3857)
  }) %>% do.call(rbind, .)

  return(regions_sf)
}

# Simulate Gaussian Process over unit grid
simulate_gp <- function(grid_res, sigma2, phi, regions, trend_means) {
  unit_grid <- unique(expand.grid(x = seq(0, 2, by = grid_res),
                           y = seq(0, 2, by = grid_res)))
  coords_mat <- as.matrix(unit_grid)
  dist_mat <- as.matrix(dist(coords_mat))
  cov_mat <- sigma2 * exp(-dist_mat / phi)
  gp_values <- MASS::mvrnorm(1, mu = rep(0, nrow(coords_mat)), Sigma = cov_mat)
  unit_grid$gp_base <- gp_values

  # Convert to sf
  gp_sf <- st_as_sf(unit_grid, coords = c("x", "y"), crs = 3857)

  # Join region polygons to the points
  gp_sf <- st_join(gp_sf, regions["name"])

  # Rename to region
  names(gp_sf)[names(gp_sf) == "name"] <- "region"

  # Add trend to base GP if region is known
  gp_sf$gp_with_trend <- gp_sf$gp_base
  has_region <- !is.na(gp_sf$region)
  gp_sf$gp_with_trend[has_region] <- gp_sf$gp_base[has_region] + trend_means[gp_sf$region[has_region]]



  list(
    gp_grid_base = gp_sf[, c("geometry", "gp_base", "region")],
    gp_grid_with_trend = gp_sf[, c("geometry", "gp_with_trend", "region")],
    coords_mat = coords_mat,
    unit_grid = unit_grid
  )
}



# Sample grid points within a region
sample_grid_points <- function(region_sf, res, n_points) {
  grid <- st_make_grid(region_sf, cellsize = res, what = "centers", square = TRUE)
  grid_sf <- st_sf(geometry = st_sfc(grid), crs = 3857)
  points_in_region <- grid_sf[region_sf, , op = st_within]

  if (nrow(points_in_region) == 0) return(NULL)

  n_sample <- min(n_points, nrow(points_in_region))
  sampled_points <- slice_sample(points_in_region, n = n_sample)
  sampled_points$region <- region_sf$name[1]
  sampled_points
}

# Simulate sampling and extract GP values
simulate_sampling <- function(regions, time_points, n_sampled, grid_res, n_grid_points, unit_grid_df, coords_mat, trend_means) {
  all_samples <- list()

  for (t in time_points) {
    sampled_names <- sample(regions$name, n_sampled)
    sampled_regions <- regions %>% filter(name %in% sampled_names)
    sampled_regions$time <- t

    for (i in seq_len(nrow(sampled_regions))) {
      region_row <- sampled_regions[i, ]
      pts <- sample_grid_points(region_row, grid_res, n_grid_points)
      if (!is.null(pts)) {
        pts$time <- t
        pts$region <- region_row$name[1]
        coords <- st_coordinates(pts)
        pts$x <- coords[, 1]
        pts$y <- coords[, 2]
        nearest_idx <- RANN::nn2(coords_mat, coords, k = 1)$nn.idx
        base_gp_values <- unit_grid_df$gp_base[nearest_idx]
        region_trend <- trend_means[pts$region]
        pts$gp_value <- base_gp_values
        pts$lp_value <- base_gp_values + region_trend
        all_samples[[length(all_samples) + 1]] <- pts
      }
    }
  }

  do.call(rbind, all_samples)
}


compute_mda_effect <- function(survey_times_data, mda_times, intervention,
                               alpha, gamma, kappa) {
  n <- length(survey_times_data)
  effect <- rep(NA, n)

  mda_effect_f <- function(v, alpha, gamma, kappa) {
    alpha*exp(-(v/gamma)^kappa)
  }
  mda_effect_f <- Vectorize(mda_effect_f, "v")

  f <- function(t, mda_times, int, alpha, gamma, kappa) {
    ind_t <- which(t > mda_times)
    u_j <- mda_times[ind_t]
    if(length(u_j) > 0) {
      out <- prod((1-mda_effect_f(t-u_j, alpha, gamma, kappa))^int[ind_t])
    } else {
      out <- 1
    }
    return(out)
  }

  f <- Vectorize(f, "t")

  for(i in 1:n) {
    effect[i] <- f(survey_times_data[i], mda_times, intervention[i,],
                   alpha, gamma, kappa)
  }
  return(effect)
}

dast_profile <- function(y, D, units_m, int_mat, survey_times_data,
                         mda_times, power_val,alpha_set) {

  p <- ncol(D)
  n <- nrow(D)

  llik_gamma <- function(par,  gamma) {
    beta <- par[1:p]
    alpha <- 1/(1+exp(-par[p+1]))

    fact <- compute_mda_effect(survey_times_data, mda_times, intervention=int_mat,
                               alpha, gamma, kappa = power_val)
    eta <- as.numeric(D%*%beta)
    prob_star <- 1/(1+exp(-eta))
    prob <- fact*prob_star

    out <- sum(y*log(prob/(1-prob)) + units_m*log(1-prob))
    return(out)
  }

  llik_alpha <- function(par,  alpha) {
    beta <- par[1:p]
    gamma <- exp(par[p+1])

    fact <- compute_mda_effect(survey_times_data, mda_times, intervention,
                               alpha, gamma, kappa = power_val)
    eta <- as.numeric(D%*%beta)
    prob_star <- 1/(1+exp(-eta))
    prob <- fact*prob_star

    out <- sum(y*log(prob/(1-prob)) + units_m*log(1-prob))
    return(out)
  }

  glm_fit <- glm(cbind(resp, den-resp) ~ -1+region, data = final_samples, family = "binomial" )
  beta_start <- coef(glm_fit)
  alphat_start <- 0
  log_gamma_start <- 0


  #alpha_set <- seq(0.8,0.999, length = acc)
  #gamma_set <- seq(1.2,1.8,length = acc)

  # alpha
  acc <- length(alpha_set)
  alpha_res <- rep(NA, acc)
  start_par <- c(beta_start, log_gamma_start)
  for(i in 1:length(alpha_set)) {
    nlminb_res <- nlminb(start_par,
                         function(x) -llik_alpha(x, alpha_set[i]))
    alpha_res[i] <- nlminb_res$objective
    start_par <- nlminb_res$par
    cat("Iter:",i,"\r")
  }

  # gamma
  #gamma_res <- rep(NA, acc)
  #start_par <- c(beta_start, alphat_start)
  #for(i in 1:length(alpha_set)) {
  #  nlminb_res <- nlminb(start_par,
  #                       function(x) -llik_gamma(x, gamma_set[i]))
  #  gamma_res[i] <- nlminb_res$objective
  #  start_par <- nlminb_res$par
  #  cat("Iter:",i,"\r")
  #}

  out <- list(
              #gamma = gamma_set,
              #val_pl_gamma = gamma_res,
              alpha = alpha_set, val_pl_alpha = alpha_res)
  return(out)

}
simulate_gp <- function(grid_res, sigma2, phi, regions = NULL, trend_means) {
  # Create a grid of unique coordinates
  unit_grid <- unique(expand.grid(x = seq(0, 2, by = grid_res),
                                  y = seq(0, 2, by = grid_res)))
  coords_mat <- as.matrix(unit_grid)

  # Construct distance and covariance matrices
  dist_mat <- as.matrix(dist(coords_mat))
  cov_mat <- sigma2 * exp(-dist_mat / phi)

  # Simulate GP values
  gp_values <- MASS::mvrnorm(1, mu = rep(0, nrow(coords_mat)), Sigma = cov_mat)
  unit_grid$gp_base <- gp_values

  # Assign regions using provided function
  unit_grid <- assign_regions_to_grid(unit_grid)

  # Add trend to GP values where region is known
  unit_grid$gp_with_trend <- unit_grid$gp_base
  has_region <- !is.na(unit_grid$region)
  unit_grid$gp_with_trend[has_region] <- unit_grid$gp_base[has_region] + trend_means[unit_grid$region[has_region]]

  # Convert to sf
  gp_sf <- st_as_sf(unit_grid, coords = c("x", "y"), crs = 3857)

  list(
    gp_grid_base = gp_sf[, c("geometry", "gp_base", "region")],
    gp_grid_with_trend = gp_sf[, c("geometry", "gp_with_trend", "region")],
    coords_mat = coords_mat,
    unit_grid = unit_grid
  )
}


find_beta_params_from_quantiles <- function(q1, p1, q2, p2, start = c(1, 1)) {
  # q1, q2 are quantile values (e.g., 0.025 and 0.975)
  # p1, p2 are the corresponding quantile probabilities (e.g., 0.1 and 0.9)

  # Objective function: sum of squared differences between target and actual quantiles
  objective <- function(par) {
    a <- par[1]
    b <- par[2]
    q1_est <- qbeta(p1, a, b)
    q2_est <- qbeta(p2, a, b)
    (q1_est - q1)^2 + (q2_est - q2)^2
  }

  # Use constrained optimization to ensure positive alpha and beta
  result <- optim(
    par = start,
    fn = objective,
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6)
  )

  list(alpha = result$par[1], beta = result$par[2])
}

fit_beta  <-  function(q1, p1, q2, p2, start_log = c(0, 0), tol = 1e-10, verbose = TRUE) {
  # Unconstrained reparameterization:
  # alpha = exp(par[1])
  # beta  = exp(par[2])

  objective <- function(par) {
    alpha <- exp(par[1])
    beta  <- exp(par[2])
    q1_est <- qbeta(p1, alpha, beta)
    q2_est <- qbeta(p2, alpha, beta)
    (q1_est - q1)^2 + (q2_est - q2)^2
  }

  result <- optim(
    par = start_log,
    fn = objective,
    method = "BFGS",  # No constraints needed now
    control = list(reltol = tol)
  )

  alpha <- exp(result$par[1])
  beta  <- exp(result$par[2])

  if (verbose) {
    cat("Estimated parameters:\n")
    cat("  alpha:", alpha, "\n")
    cat("  beta: ", beta, "\n")
    cat("  qbeta(", p1, ") = ", qbeta(p1, alpha, beta), "\n")
    cat("  qbeta(", p2, ") = ", qbeta(p2, alpha, beta), "\n")
  }

  return(list(alpha = alpha, beta = beta, value = result$value, convergence = result$convergence))
}

assign_regions_to_grid <- function(grid) {
  grid <- as.data.frame(grid)  # Ensure it's a data frame
  grid$region <- with(grid, ifelse(
    x >= 0 & x <= 1 & y >= 0 & y <= 1, "A",
    ifelse(
      x > 1 & x <= 2 & y >= 0 & y <= 1, "B",
      ifelse(
        x >= 0 & x <= 1 & y > 1 & y <= 2, "C",
        ifelse(
          x > 1 & x <= 2 & y > 1 & y <= 2, "D",
          NA
        )
      )
    )
  ))
  return(grid)
}

