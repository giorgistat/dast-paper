#!/usr/bin/env Rscript
# =============================================================================
# STH in Kenya — DAST application (tidied for GitHub)
# =============================================================================
# Purpose
#   Fit the DAST model for three STH species (Ascaris, Trichuris, Hookworm)
#   using Kenyan survey data and subcounty-level MDA coverage.
#   Produce MDA effect comparison plot and animated prevalence maps (grid and
#   subcounty aggregations).
#
# Inputs (expected in working directory)
#   - KenyaData.csv                      (survey data with x, y, year, subcounty, noStudSamp, noASC, noTT, noHKW)
#   - KenyaCoverage_TotalMDA.shp (+ sidecar files)  (polygon shapefile with 'subcounty' and covYYYY fields)
#
# Outputs (written to ./outputs/)
#   - MDA_effects.png
#   - STH_grid_animation.gif
#   - STH_admin_animation.gif
#
# Notes
#   - Requires RiskMap >= 1.0.0 with functions: dast(), pred_over_grid(),
#     pred_target_grid(), pred_target_shp(), plot_mda(), create_grid().
#   - Computation can be heavy; adjust image resolution, ranges, and years if needed.
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

# --------------------------- 0. Libraries -----------------------------------
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(RiskMap)
  library(terra)
  library(magick)
  library(cowplot)
})

if (!dir.exists("outputs")) dir.create("outputs")

# --------------------------- 1. Data ----------------------------------------
# 1.1 Survey data (point locations in ESRI:102022 as provided, then projected to UTM for Kenya)
sth <- read.csv("KenyaData.csv")

stopifnot(all(c("x","y","year","subcounty","noStudSamp","noASC","noTT","noHKW") %in% names(sth)))

sth <- st_as_sf(sth, coords = c("x", "y"), crs = "ESRI:102022")
data_utm <- propose_utm(sth)           # choose a suitable UTM CRS for Kenya
sth <- st_transform(sth, crs = data_utm)

# 1.2 Subcounty MDA polygons
mda_data_sf <- st_read("KenyaCoverage_TotalMDA.shp", quiet = TRUE)
mda_data_sf <- st_transform(mda_data_sf, crs = data_utm)
mda_data <- as.data.frame(mda_data_sf)

stopifnot("subcounty" %in% names(mda_data))

# --------------------------- 2. Preprocess MDA ------------------------------
# 2.1 Identify MDA coverage columns (covYYYY) and extract years
mda_columns <- grepl("^cov", names(mda_data), ignore.case = TRUE)
cov_columns <- names(mda_data)[mda_columns]

if (length(cov_columns) == 0) stop("No 'covYYYY' columns found in shapefile attributes.")

# Years are given as suffix after 'cov' and represent years since 2000
mda_times <- sort(as.numeric(sub("cov", "", tolower(cov_columns)))) + 2000

# 2.2 Build the intervention matrix (n_locations × n_mda_years)
n_loc <- nrow(sth)
n_mda <- length(mda_times)
intervention <- matrix(0L, nrow = n_loc, ncol = n_mda)

# Ensure a safe character match on 'subcounty'
sth$subcounty <- as.character(sth$subcounty)
mda_data$subcounty <- as.character(mda_data$subcounty)

for (i in seq_len(n_loc)) {
  idx <- which(mda_data$subcounty == sth$subcounty[i])
  if (length(idx) == 1) {
    intervention[i, ] <- as.integer(mda_data[idx, cov_columns] > 0)
  } else {
    # If no unique match, leave zeros (assume no MDA info)
    intervention[i, ] <- 0L
  }
}

# 2.3 Add planar coordinates explicitly for gp() terms
sth$x <- st_coordinates(sth)[, 1]
sth$y <- st_coordinates(sth)[, 2]

# --------------------------- 3. Fit DAST per species ------------------------
# Parameter starting values (example values taken from prior fits)
par0_asc <- list(beta = -2.44252,
                 sigma2 = exp(3.06358),
                 phi = exp(5.25905),
                 alpha = 1/(1 + exp(-(-0.878661))),
                 gamma = exp(2.80133),
                 tau2 = 1)

par0_tt  <- list(beta = -2.40621,
                 sigma2 = exp(3.58584),
                 phi = exp(5.37115),
                 alpha = 1/(1 + exp(-(-0.626555))),
                 gamma = exp(1.71385))

par0_hkw <- list(beta = -1.98672,
                 sigma2 = exp(3.45462),
                 phi = exp(5.43809),
                 alpha = 1/(1 + exp(-(2.53378))),
                 gamma = exp(1.65223))

# 3.1 Ascaris: spatio-temporal GP
dast_asc <-
  dast(formula = noASC ~ gp(x, y, year, nugget = NULL),
       den = noStudSamp, power_val = 1,
       data = sth,
       par0 = par0_asc, start_pars = par0_asc,
       mda_times = mda_times,
       int_mat = intervention)

# 3.2 Trichuris: simpler GP (survey_times provided separately)
dast_tt <-
  dast(formula = noTT ~ gp(),
       den = noStudSamp, power_val = 1,
       data = sth,
       par0 = par0_tt, start_pars = par0_tt,
       survey_times = year,
       mda_times = mda_times,
       int_mat = intervention)

# 3.3 Hookworm: simpler GP
dast_hkw <-
  dast(formula = noHKW ~ gp(),
       den = noStudSamp, power_val = 1,
       data = sth,
       par0 = par0_hkw, start_pars = par0_hkw,
       survey_times = year,
       mda_times = mda_times,
       int_mat = intervention)

# Print concise summaries
print(summary(dast_asc))
print(summary(dast_tt))
print(summary(dast_hkw))

# --------------------------- 4. MDA effect plot -----------------------------
# Compare the estimated MDA impact functions across species.
gg_asc <- plot_mda(dast_asc, upper_f = 1, lower_f = 0, conf_level = 0.95, x_max = 20)
gg_tt  <- plot_mda(dast_tt,  upper_f = 1, lower_f = 0, conf_level = 0.95, x_max = 20)
gg_hkw <- plot_mda(dast_hkw, upper_f = 1, lower_f = 0, conf_level = 0.95, x_max = 20)

# Extract plotted data from each ggplot object and label them
df_asc <- ggplot_build(gg_asc)$data[[1]] %>% mutate(Species = "Ascaris")
df_tt  <- ggplot_build(gg_tt)$data[[1]]  %>% mutate(Species = "Trichuris")
df_hkw <- ggplot_build(gg_hkw)$data[[1]] %>% mutate(Species = "Hookworm")

df_all <- dplyr::bind_rows(df_asc, df_tt, df_hkw)

gg_combined <- ggplot(df_all, aes(x = x, y = y, colour = Species)) +
  geom_line(linewidth = 1) +
  theme_minimal(base_size = 13) +
  labs(x = "Years since MDA", y = "Relative prevalence reduction",
       title = "Estimated MDA impact functions (Kenya, STH)") +
  scale_colour_manual(values = c("Ascaris" = "blue", "Trichuris" = "red", "Hookworm" = "darkgreen")) +
  theme(legend.position = "bottom")

ggsave("outputs/MDA_effects.png", gg_combined, width = 7, height = 5, dpi = 300)

# --------------------------- 5. Predictions over grid -----------------------
# 5.1 Create a regular prediction grid over the study area
grid_pred <- create_grid(st_as_sf(st_union(mda_data_sf)), spat_res = 5)

# 5.2 Predict over grid (joint summaries)
pred_dast_asc <- pred_over_grid(dast_asc, grid_pred = grid_pred, type = "joint")
pred_dast_tt  <- pred_over_grid(dast_tt,  grid_pred = grid_pred, type = "joint")
pred_dast_hkw <- pred_over_grid(dast_hkw, grid_pred = grid_pred, type = "joint")

# --------------------------- 6. Attach MDA to grid --------------------------
# Join polygon attributes to grid cells to get MDA histories per cell
grid_pred_ext <- st_join(st_as_sf(pred_dast_asc$grid_pred), mda_data_sf, join = st_within)
mda_grid <- st_drop_geometry(grid_pred_ext[, cov_columns, drop = FALSE])

# --------------------------- 7. Time-varying predictions --------------------
# Define prediction years
time_pred <- seq(2012, 2023, by = 1)
n_time <- length(time_pred)

pred_dast_grid_asc <- vector("list", n_time)
pred_dast_grid_tt  <- vector("list", n_time)
pred_dast_grid_hkw <- vector("list", n_time)

pred_dast_grid_asc_sc <- vector("list", n_time)
pred_dast_grid_tt_sc  <- vector("list", n_time)
pred_dast_grid_hkw_sc <- vector("list", n_time)

for (i in seq_len(n_time)) {
  # Grid-scale predictions
  pred_dast_grid_asc[[i]] <-
    pred_target_grid(dast_asc, mda_grid = mda_grid, time_pred = time_pred[i],
                     f_target = list(prev = function(x) 1/(1+exp(-x))),
                     pd_summary = list(mean = mean))

  pred_dast_grid_tt[[i]] <-
    pred_target_grid(dast_tt,  mda_grid = mda_grid, time_pred = time_pred[i],
                     f_target = list(prev = function(x) 1/(1+exp(-x))),
                     pd_summary = list(mean = mean))

  pred_dast_grid_hkw[[i]] <-
    pred_target_grid(dast_hkw, mda_grid = mda_grid, time_pred = time_pred[i],
                     f_target = list(prev = function(x) 1/(1+exp(-x))),
                     pd_summary = list(mean = mean))

  # Subcounty aggregation (shapefile level)
  pred_dast_grid_asc_sc[[i]] <-
    pred_target_shp(dast_asc, shp = mda_data_sf, col_names = "subcounty",
                    mda_grid = mda_grid, time_pred = time_pred[i],
                    shp_target = mean,
                    f_target = list(prev = function(x) 1/(1+exp(-x))),
                    pd_summary = list(mean = mean))

  pred_dast_grid_tt_sc[[i]] <-
    pred_target_shp(dast_tt, shp = mda_data_sf, col_names = "subcounty",
                    mda_grid = mda_grid, time_pred = time_pred[i],
                    shp_target = mean,
                    f_target = list(prev = function(x) 1/(1+exp(-x))),
                    pd_summary = list(mean = mean))

  pred_dast_grid_hkw_sc[[i]] <-
    pred_target_shp(dast_hkw, shp = mda_data_sf, col_names = "subcounty",
                    mda_grid = mda_grid, time_pred = time_pred[i],
                    shp_target = mean,
                    f_target = list(prev = function(x) 1/(1+exp(-x))),
                    pd_summary = list(mean = mean))
}

# --------------------------- 8. GIF: grid-level maps ------------------------
dir.create("gif_frames", showWarnings = FALSE)

for (i in seq_len(n_time)) {
  png_filename <- sprintf("gif_frames/frame_%03d.png", i)
  png(png_filename, width = 1800, height = 600, res = 150)
  par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))

  plot(pred_dast_grid_asc[[i]], which_target = "prev", which_summary = "mean",
       main = paste("Ascaris: Year", time_pred[i]), range = c(0, 0.5), cex.main = 2)

  plot(pred_dast_grid_tt[[i]], which_target = "prev", which_summary = "mean",
       main = paste("Trichuris: Year", time_pred[i]), range = c(0, 0.85), cex.main = 2)

  plot(pred_dast_grid_hkw[[i]], which_target = "prev", which_summary = "mean",
       main = paste("Hookworm: Year", time_pred[i]), range = c(0, 0.4), cex.main = 2)

  dev.off()
}

frames <- list.files("gif_frames", pattern = "png$", full.names = TRUE)
animation <- image_read(frames)
animation <- image_animate(animation, fps = 2)
image_write(animation, "outputs/STH_grid_animation.gif")

# --------------------------- 9. GIF: subcounty maps -------------------------
dir.create("gif_frames_sc", showWarnings = FALSE)

for (i in seq_len(n_time)) {
  year_label <- time_pred[i]

  p1 <- plot(pred_dast_grid_asc_sc[[i]], which_target = "prev", which_summary = "mean") +
    ggtitle(paste("Ascaris: Year", year_label)) +
    scale_fill_viridis_c(limits = c(0, 0.35), na.value = "white", name = "Prev") +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

  p2 <- plot(pred_dast_grid_tt_sc[[i]], which_target = "prev", which_summary = "mean") +
    ggtitle(paste("Trichuris: Year", year_label)) +
    scale_fill_viridis_c(limits = c(0, 0.25), na.value = "white", name = "Prev") +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

  p3 <- plot(pred_dast_grid_hkw_sc[[i]], which_target = "prev", which_summary = "mean") +
    ggtitle(paste("Hookworm: Year", year_label)) +
    scale_fill_viridis_c(limits = c(0, 0.25), na.value = "white", name = "Prev") +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

  combined <- cowplot::plot_grid(p1, p2, p3, ncol = 3, align = "hv", rel_widths = c(1, 1, 1))

  png_filename <- sprintf("gif_frames_sc/frame_%03d.png", i)
  ggsave(png_filename, combined, width = 18, height = 6, dpi = 150)
}

frames2 <- list.files("gif_frames_sc", pattern = "png$", full.names = TRUE)
animation2 <- image_read(frames2)
animation2 <- image_animate(animation2, fps = 2)
image_write(animation2, "outputs/STH_admin_animation.gif")

# --------------------------- 10. Optional: AnPIT figure & CRPS table --------
# The following section assumes that the objects u_val, m_AnPIT_* and CRPS_* are
# present in the environment from prior computations. It will be skipped otherwise.

if (exists("u_val") &&
    all(c("m_AnPIT_asc_st","m_AnPIT_asc_s","m_AnPIT_asc_dast",
          "m_AnPIT_tt_st","m_AnPIT_tt_s","m_AnPIT_tt_dast",
          "m_AnPIT_hkw_st","m_AnPIT_hkw_s","m_AnPIT_hkw_dast") %in% ls())) {

  suppressPackageStartupMessages({
    library(tidyr); library(patchwork); library(knitr)
  })

  df_plot <- tibble::tibble(
    u        = rep(u_val, times = 9),
    AnPIT    = c(m_AnPIT_asc_st,  m_AnPIT_asc_s,  m_AnPIT_asc_dast,
                 m_AnPIT_tt_st,   m_AnPIT_tt_s,   m_AnPIT_tt_dast,
                 m_AnPIT_hkw_st,  m_AnPIT_hkw_s,  m_AnPIT_hkw_dast),
    Model    = rep(rep(c("ST", "S", "DAST"), each = length(u_val)), times = 3),
    Parasite = rep(c("Ascaris", "Trichuris", "Hookworm"),
                   each = length(u_val) * 3)
  )

  anpit_plot <- ggplot(df_plot, aes(x = u, y = AnPIT, colour = Model)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, colour = "black") +
    facet_wrap(~Parasite, ncol = 1, strip.position = "left") +
    scale_colour_manual(values = c(ST = "blue", S = "forestgreen", DAST = "red")) +
    labs(x = "u", y = "Average non-randomized PIT (AnPIT)", colour = NULL) +
    theme_minimal() + theme(legend.position = "bottom", strip.placement = "outside")

  ggsave("outputs/AnPIT_sth.pdf", anpit_plot, width = 6, height = 10)

  if (all(c("CRPS_asc_st","CRPS_tt_st","CRPS_hkw_st",
            "CRPS_asc_s","CRPS_tt_s","CRPS_hkw_s",
            "CRPS_asc_dast","CRPS_tt_dast","CRPS_hkw_dast") %in% ls())) {

    mean_crps <- matrix(
      c(mean(CRPS_asc_st,  na.rm = TRUE), mean(CRPS_tt_st,  na.rm = TRUE), mean(CRPS_hkw_st,  na.rm = TRUE),
        mean(CRPS_asc_s,   na.rm = TRUE), mean(CRPS_tt_s,   na.rm = TRUE), mean(CRPS_hkw_s,   na.rm = TRUE),
        mean(CRPS_asc_dast,na.rm = TRUE), mean(CRPS_tt_dast,na.rm = TRUE), mean(CRPS_hkw_dast,na.rm = TRUE)),
      nrow = 3, byrow = TRUE,
      dimnames = list(Model     = c("ST", "S", "DAST"),
                      Parasite  = c("Ascaris", "Trichuris", "Hookworm"))
    )

    print(knitr::kable(round(mean_crps, 3),
                       format   = "latex",
                       booktabs = TRUE,
                       caption  = "Mean CRPS on the held-out data (lower is better)"))
  }
}

message("STH Kenya DAST application finished. See 'outputs/' for figures and GIFs.")
