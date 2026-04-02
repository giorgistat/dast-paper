# =============================================================================
# Figures produced in the paper
# =============================================================================
# Code used to produce some the figures reported in the paper 
#
# Inputs (expected in ./data/)
#   - KenyaData.csv (survey data for STH Kenya)
#   - KenyaCoverage_TotalMDA.shp (subcounty shapefile with MDA data)
#   - lf_madagascar.csv (survey data for LF Madagascar)
#   - lf_MDA_madagascar.csv (MDA data for LF)
#   - mdg_adm0.shp and mdg_adm2.shp (shapefiles for Madagascar)
#   - simulations.qs (results from the simulation study)
#
# Outputs (written to ./figures/)
#   - figure1_schematic.pdf               (Fig. 1)
#   - sim_example.pdf                     (Fig. 2)
#   - mda_impact.pdf                      (Fig. 3)
#   - sim_bias.pdf                        (Fig. 4)
#   - sim_RMSE.pdf                        (Fig. 5)
#   - sth_kenya_panels.pdf                (Fig. 6)
#   - combined_MDA_survey_timeline.pdf    (Fig. 7)
#   - cum_mda_kenya.pdf                   (Fig. 8)
#   - mda_impact_sth_with_CI.pdf          (Fig. 9)
#   - AnPIT_sth.pdf                       (Fig. 10)
#   - lf_madagascar_panels.pdf            (Fig. 11)
#   - combined_MDA_survey_timeline_lf.pdf (Fig. 12)
#   - cum_mda_mdg.pdf                     (Fig. 13)
#   - lf_anpit.pdf                        (Fig. 14)
#   - lf_rounds.pdf                       (Fig. 15)

# =============================================================================#

# 0. LIBRARIES -----------------------------------------------------------------

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(scales)
  library(ggspatial)
  library(rgeoboundaries)
})

if (!dir.exists("outputs")) dir.create("outputs")

# FIGURE 1. --------------------------------------------------------------------

## STEP 1: SIMULATE REGIONS AND SUBREGIONS -------------------------------------

# Define coordinates of sub-regions (squares)
district_coords <- list(
  A = matrix(c(0, 1, 1, 0, 0, 1, 1, 2, 2, 1), ncol = 2),
  B = matrix(c(1, 2, 2, 1, 1, 1, 1, 2, 2, 1), ncol = 2),
  C = matrix(c(1, 2, 2, 1, 1, 0, 0, 1, 1, 0), ncol = 2),
  D = matrix(c(0, 1, 1, 0, 0, 0, 0, 1, 1, 0), ncol = 2)
)

# Create sf polygons
districts_sf <- lapply(names(district_coords), function(name) {
  st_polygon(list(district_coords[[name]])) %>%
    st_sfc(crs = 4326) %>%
    st_sf(geometry = ., district = name)
}) %>%
  do.call(rbind, .)

# Plot the region
ggplot(data = districts_sf) +
  geom_sf(aes(fill = district)) +
  geom_sf_label(data = districts_sf, aes(label = district)) 

## STEP 2: SIMULATE SPTIAL GP --------------------------------------------------

# Sample n points in the region
n_points <- 200
set.seed(123)
sample_pts <- districts_sf |> 
  st_sample(size = n_points, type = "random", exact = TRUE) |> 
  st_as_sf()

# Extract coordinates for covariance calculation
coords <- st_coordinates(sample_pts)

# Spatial GP parameters
sigma2 <- 1     # Variance
phi <- 0.5      # Range

# Distance matrix
D <- as.matrix(dist(coords))

# Exponential covariance matrix
Sigma <- sigma2 * exp(-D / phi)

# Simulate spatial effect
set.seed(123)
S <- MASS::mvrnorm(1, mu = rep(0, n_points), Sigma = Sigma)

# Attach S(x) to points
sample_pts$S <- S

# Assign districts using spatial join
sample_pts <- sample_pts |> 
  st_join(districts_sf)

# Visualize simulated spatial field
ggplot(sample_pts) +
  geom_sf(aes(color = S)) +
  geom_sf(data = districts_sf, fill = NA, col = "black") +
  geom_sf_label(data = districts_sf, aes(label = district)) +
  scale_color_gradient2(low = "blue", high = "red") +
  theme_minimal() +
  ggtitle("Simulated Spatial Gaussian Process at Random Locations")

## STEP 3: ADD DISTRICT EFFECTS ------------------------------------------------

# Define district effects (increasing from A to D)
district_effects <- c(A = -1, B = 0, C = 1, D = 3.5)

# Assign district effect
sample_pts$beta <- district_effects[sample_pts$district]

# Compute counterfactual prevalence using inverse logit
inv_logit <- function(x) exp(x) / (1 + exp(x))
sample_pts$P_star <- inv_logit(sample_pts$beta + sample_pts$S)

sample_pts |> 
  group_by(district) |> 
  summarise(mean(P_star))

# Plot counterfactual prevalence
ggplot(sample_pts) +
  geom_sf(aes(color = P_star)) +
  geom_sf_label(data = districts_sf, aes(label = district)) +
  geom_sf(data = districts_sf, fill = NA, col = "black") +
  viridis::scale_color_viridis(option = "C", name = "P*(x)") +
  theme_minimal() +
  ggtitle("Counterfactual Prevalence P*(x) by District and Spatial Field")

## STEP 4: APPLY MDA EFFECT OVER TIME ------------------------------------------

# Time points: t = 0 (baseline), ..., t = 4 (before 5th MDA)
t <- 3
time_points <- 0:t

# MDA decay parameters
alpha <- 0.25
gamma <- 1.5
kappa <- 0.5

# MDA effect function
f_mda <- function(v, alpha, gamma, kappa) {
  alpha * exp(- (v / gamma)^kappa)
}

# Compute f(v) for v = 0 to 3 (for t = 0, 1, 2, 3, 4)
lags <- 0:(t - 1)
f_vals <- f_mda(lags, alpha = alpha, gamma = gamma, kappa = kappa)

# Plot
mda_df <- data.frame(lag = lags, f = f_vals)
ggplot(mda_df, aes(x = lag, y = f)) +
  geom_line(linewidth = 1.2, color = "darkred") +
  geom_point(size = 3, color = "darkred") +
  scale_x_continuous(breaks = lags) +
  labs(title = "MDA Effect Function f(v)",
       x = "Lag (v = t - j)",
       y = "MDA Effect (f(v))") +
  theme_minimal()

# MDA multiplier at each t (surveys conducted BEFORE MDA rounds)
mda_decay <- sapply(time_points, function(t) {
  if (t == 0) {
    1  # No MDA effect yet
  } else {
    prod(1 - sapply(0:(t - 1), function(v) f_mda(v, alpha, gamma, kappa)))
  }
})
names(mda_decay) <- paste0("t", time_points)

# Expand sample points across all time points
sample_pts[c("lon", "lat")] <- st_coordinates(sample_pts)
pts_long <- sample_pts |> 
  st_drop_geometry() |> 
  mutate(id = row_number()) |> 
  tidyr::crossing(time = time_points)

# Assign prevalence with MDA effect
pts_long$P <- pts_long$P_star * as.numeric(mda_decay[paste0("t", pts_long$time)])

# Check district prevalence over time
avg_prev_district <- pts_long |> 
  group_by(time, district) |> 
  summarise(mean_prev = mean(P))

avg_prev_district

districts_time <- inner_join(districts_sf, avg_prev_district) 

# Create sampled data
# Train set
pts_long_train <- pts_long |> 
  filter((district == "A" & time == 0) |
           (district == "B" & time == 1) |
           (district == "C" & time == 2) |
           (district == "D" & time == 3)
  ) |> 
  group_by(district) |> 
  sample_n(10)

# Custom labels for time points and districts
districts_time$time_label <- paste0("bold(t[", districts_time$time, "])")
pts_long_train$time_label <- paste0("bold(t[", pts_long_train$time, "])")

# Define the levels and their corresponding subscript labels
district_levels <- c("A", "B", "C", "D")
district_labels <- paste0("U[", seq_along(district_levels), "]")

# Create a factor with custom labels (subscripts as expressions)
districts_time$district_label <- factor(
  districts_time$district,
  levels = district_levels,
  labels = district_labels
)

g1 <- ggplot(pts_long_train) +
  geom_sf(data = districts_time, aes(fill = mean_prev), col = "black") +
  geom_point(aes(x = lon, y = lat)) +
  geom_sf_label(data = districts_time, aes(label = district_label), parse = T, fill = "white") +
  scale_fill_distiller("Prevalence", palette = 18, direction = 1, 
                       breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  facet_wrap(~ time_label, ncol = 4, labeller = label_parsed) +
  theme_void() +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  theme(legend.position = "top", legend.key.width = unit(1.5, "cm"))
g1

# Regional average true
avg_prev_region_true <- pts_long |> 
  group_by(time) |> 
  summarise(mean_prev = mean(P), .groups = "drop") |> 
  mutate(type = "True")

# Regional average biased
avg_prev_region_bias <- pts_long_train |> 
  group_by(time) |> 
  summarise(mean_prev = mean(P), .groups = "drop") |> 
  mutate(type = "Empirical")

# Combine together and plot
avg_prev_region <- rbind(avg_prev_region_bias, avg_prev_region_true)

# Get the last time point for each type to position the label
label_data <- avg_prev_region  |> 
  group_by(type) |> 
  filter(time == max(time))

g2 <- ggplot(avg_prev_region, aes(x = time, y = mean_prev, group = type)) +
  geom_line(linewidth = 1.2, aes(linetype = type)) +
  geom_point(data = avg_prev_region_bias, size = 2) +
  geom_text(data = label_data,
            aes(label = type),
            hjust = 0.5, vjust = -0.7, size = 4.5) +
  labs(x = "Time (years)", y = "Mean Prevalence",
       title = "Regional Average Prevalence Over Time") +
  scale_linetype(NULL) +
  scale_x_continuous(labels = function(x) parse(text = paste0("bold(t[", x, "])")),
                     limits = c(0, max(avg_prev_region$time) + 0.08)) +
  theme_minimal() +
  theme(legend.position = "none") 
g2

# Combine plots vertically with labels
cowplot::plot_grid(g1, g2, ncol = 1,labels = c("A", "B"), 
                   label_size = 14, align = "v", rel_heights = c(0.45, 0.55))
ggsave("figures/figure1_schematic.pdf", width = 10, height = 7)

# FIGURE 2. --------------------------------------------------------------------


# FIGURE 3. --------------------------------------------------------------------

sim_df <- qs::qread("simulations.qs")

sim_df <- sim_df |>
  dplyr::mutate(
    scenario = factor(scenario, levels = 1:2,
                      labels = c("Scenario 1", "Scenario 2")),
    penalty  = factor(penalty, levels = c("no", "yes"),
                      labels = c("Unpenalised", "Penalised")),
    time_f   = factor(time, levels = paste0("T", 1:5),
                      labels = paste0("t=", 0:4))
  ) |>
  dplyr::filter(!is.na(pred_prev))

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Shared intermediate computations
# ─────────────────────────────────────────────────────────────────────────────

# sd of predicted logits per (region, scenario, penalty, time)
cell_scale <- sim_df |>
  dplyr::group_by(region, scenario, penalty, time) |>
  dplyr::summarise(S_L = sd(pred_logit, na.rm = TRUE), .groups = "drop")

# Per-observation standardised error
std_obs <- sim_df |>
  dplyr::mutate(deriv_local = pred_prev * (1 - pred_prev)) |>
  dplyr::inner_join(cell_scale, by = c("region", "scenario", "penalty", "time")) |>
  dplyr::mutate(zP = (pred_prev - true_prev) / (deriv_local * S_L))

# Average over regions within each replicate
std_rep <- std_obs |>
  dplyr::group_by(scenario, penalty, time_f, sim) |>
  dplyr::summarise(
    zP     = mean(zP,    na.rm = TRUE),
    zP_ms  = mean(zP^2,  na.rm = TRUE),
    zP_ma  = mean(abs(zP), na.rm = TRUE),
    .groups = "drop"
  )

# Summarise across replicates
std_sum <- std_rep |>
  dplyr::group_by(scenario, penalty, time_f) |>
  dplyr::summarise(
    SBias     = mean(zP,    na.rm = TRUE),
    SBias_q25 = quantile(zP, 0.25, na.rm = TRUE),
    SBias_q75 = quantile(zP, 0.75, na.rm = TRUE),
    sRMSE     = sqrt(mean(zP_ms, na.rm = TRUE)),
    .groups   = "drop"
  )

# IU-level mean e_{k,t} across replicates
mean_e_iu <- std_obs |>
  dplyr::group_by(region, scenario, penalty, time_f) |>
  dplyr::summarise(
    mean_e = mean(zP, na.rm = TRUE),
    .groups = "drop"
  )

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Shared theme and palette
# ─────────────────────────────────────────────────────────────────────────────

pal      <- c("Unpenalised" = "#5DA5DA", "Penalised" = "#F15854")
t4_col   <- which(levels(factor(paste0("t=", 0:4))) == "t=4")

theme_pub <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "top"
  )

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Figure: sBias over time
# ─────────────────────────────────────────────────────────────────────────────

sBias <- ggplot(std_sum,
                aes(x = time_f, y = SBias, colour = penalty, group = penalty)) +
  annotate("rect",
           xmin = t4_col - 0.5, xmax = t4_col + 0.5,
           ymin = -Inf, ymax = Inf, alpha = 0.08) +
  geom_errorbar(aes(ymin = SBias_q25, ymax = SBias_q75),
                width = 0.15, linewidth = 0.4, alpha = 0.6) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8, shape = 21, fill = "white", stroke = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~scenario, nrow = 1) +
  scale_colour_manual("", values = pal) +
  labs(x = NULL, y = "sBias") +
  theme_pub

ggsave("figures/sim_bias.pdf", sBias, width = 9, height = 5, dpi = 300)
message("Saved figures/sim_bias.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# 5.  Figure: sRMSE over time
# ─────────────────────────────────────────────────────────────────────────────

sRMSE <- ggplot(std_sum,
                aes(x = time_f, y = sRMSE, colour = penalty, group = penalty)) +
  annotate("rect",
           xmin = t4_col - 0.5, xmax = t4_col + 0.5,
           ymin = -Inf, ymax = Inf, alpha = 0.08) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  facet_wrap(~scenario, nrow = 1) +
  scale_colour_manual("", values = pal) +
  labs(x = NULL, y = "sRMSE") +
  theme_pub

ggsave("figures/sim_RMSE.pdf", sRMSE, width = 9, height = 5, dpi = 300)
message("Saved figures/sim_RMSE.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# 6.  Figure: IU-level choropleth of mean e_{k,t}
#
#   The four IUs (A, B, C, D) are arranged in a 2x2 grid on [0,2]x[0,2].
#   Adjust the polygon coordinates below if your simulation uses a different
#   region layout.
#     A = [0,1]x[0,1]   B = [1,2]x[0,1]
#     C = [0,1]x[1,2]   D = [1,2]x[1,2]
# ─────────────────────────────────────────────────────────────────────────────

iu_boxes <- list(
  R1 = sf::st_polygon(list(rbind(c(0,0),c(1,0),c(1,1),c(0,1),c(0,0)))),
  R2 = sf::st_polygon(list(rbind(c(1,0),c(2,0),c(2,1),c(1,1),c(1,0)))),
  R3 = sf::st_polygon(list(rbind(c(0,1),c(1,1),c(1,2),c(0,2),c(0,1)))),
  R4 = sf::st_polygon(list(rbind(c(1,1),c(2,1),c(2,2),c(1,2),c(1,1))))
)

iu_sf <- sf::st_sf(
  region   = names(iu_boxes),
  geometry = sf::st_sfc(iu_boxes),
  crs      = NA_crs_
) |>
  dplyr::mutate(
    x_cent = c(0.5, 1.5, 0.5, 1.5),
    y_cent = c(0.5, 0.5, 1.5, 1.5)
  )

plot_sf <- iu_sf |>
  dplyr::inner_join(mean_e_iu, by = "region")

lim <- range(plot_sf$mean_e)

choro <- ggplot(plot_sf) +
  geom_sf(aes(fill = mean_e), colour = "white", linewidth = 0.7) +
  geom_text(aes(x = x_cent, y = y_cent, label = region),
            size = 3.5, fontface = "bold", colour = "grey20") +
  facet_grid(rows = vars(scenario, penalty),
             cols = vars(time_f)) +
  scale_fill_gradient2(
    low      = "#d6604d",
    mid      = "white",
    high     = "#4393c3",
    midpoint = sum(lim)/2,
    limits   = lim,
    name     = expression(bar(e)[kt])
  ) +
  labs(
    title   = expression("Mean standardised prediction error by IU and time point"),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    strip.text      = element_text(face = "bold", size = 9),
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 12),
    plot.caption    = element_text(size = 9, colour = "grey40")
  )

ggsave("figures/sim_choropleth_e.pdf", choro,
       width = 14, height = 9, dpi = 300)
message("Saved figures/sim_choropleth_e.pdf")

# FIGURE 6. --------------------------------------------------------------------

## 1. IMPORT DATA --------------------------------------------------------------

# Survey data STH (point locations in ESRI:102022 as provided)
sth <- read.csv("data/KenyaData.csv")

# Kenya boudaries
kenya_bg <- gb_adm0("Kenya")

# Kenya subcounty MDA polygons
subcounty <- st_read("data/KenyaCoverage_TotalMDA.shp", quiet = TRUE)

## 2. DATA PROCESSING ----------------------------------------------------------

# Set and change CRS
sth <- st_as_sf(sth, coords = c("x", "y"), crs = "ESRI:102022")
sth <- st_transform(sth, crs = 4326)
subcounty <- st_transform(subcounty, crs = 4326)

# Define time panels STH
panel_labs_sth <- c("2012 (Baseline)", "2013–2014", "2015–2016",
                    "2017–2018", "2021–2022")

sth_pan <- sth |>
  mutate(
    panel = case_when(
      year == 2012 ~ "2012 (Baseline)",
      year %in% 2013:2014 ~ "2013–2014",
      year %in% 2015:2016 ~ "2015–2016",
      year %in% 2017:2018 ~ "2017–2018",
      year %in% 2021:2022 ~ "2021–2022",
      .default = NA_character_
    ),
    panel = factor(panel, levels = panel_labs_sth)
  )

## 3. PREVALENCE MAP -----------------------------------------------------------

base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "top"
  )

# A function to plot a single year-range panel
plot_panel <- function(lbl) {
  pdat <- filter(sth_pan, panel == lbl)
  ggplot() +
    geom_sf(data = kenya_bg, fill = "grey95", color = "black") +
    geom_sf(
      data = pdat,
      aes(fill = prevSTH),
      stroke = 0.2, shape = 21, size = 1.1, alpha = .95
    ) +
    scale_fill_distiller(
      type = "div",
      palette = "RdYlBu",
      direction = -1,
      limits = c(0, 0.75),
      breaks = seq(0, 0.75, by = 0.1),
      labels = scales::label_percent(accuracy = 1),
      name = "Prevalence (%)",
      guide = guide_colorbar(
        title.position = "top",        
        title.hjust = 0.5,             
        barwidth = unit(6, "in"),      
        barheight = unit(0.4, "cm")
      )) +
    coord_sf() +
    labs(title = lbl) +
    base_theme
}

# Build each panel
p_2012 <- plot_panel(panel_labs_sth[1])
p_1314 <- plot_panel(panel_labs_sth[2])
p_1516 <- plot_panel(panel_labs_sth[3])
p_1718 <- plot_panel(panel_labs_sth[4])
p_2122 <- plot_panel(panel_labs_sth[5])

# Subcounty-only panel
p_subc <- ggplot() +
  geom_sf(data = kenya_bg, fill = "grey95", color = NA) +
  geom_sf(data = subcounty, fill = "grey", color = "grey40", linewidth = 0.15) +
  labs(title = "Subcounties") +
  base_theme +
  theme(legend.position = "none") +
  annotation_scale(location = "bl", width_hint = 0.25) +  
  annotation_north_arrow(location = "tr", which_north = "true", 
                         height = unit(.5, "cm"), width = unit(.5, "cm")) 

# Assemble 3×2 layout and add title/caption
plt <- (guide_area() / (
  p_2012 | p_1314 | p_1516
) /
  (p_1718 | p_2122 | p_subc)) +
  plot_layout(guides = "collect", heights = unit(c(1, 1, 1), c("cm", "null", "null"))) &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title = element_text(size = 11, face = "bold", vjust = 0.3),
    legend.text = element_text(size = 9),
  )

plt
ggsave(filename = "figures/sth_kenya_panels.pdf",
       plot = plt, device = cairo_pdf, 
       bg = "white", width = 8, height = 7.5)


# FIGURE 7. --------------------------------------------------------------------

# FIGURE 8. --------------------------------------------------------------------

mda_data <- st_drop_geometry(subcounty)

#  Identify MDA coverage columns (covYYYY) and extract years
mda_columns <- grepl("^cov", names(mda_data), ignore.case = TRUE)
cov_columns <- names(mda_data)[mda_columns]

# Years are given as suffix after 'cov' and represent years since 2000
mda_times <- sort(as.numeric(sub("cov", "", tolower(cov_columns)))) + 2000
names(mda_data)[mda_columns] <- mda_times

# Compute cum  mda
mda_cum <- mda_data |>
  as_tibble() |> 
  select(subcounty, `2013`:`2023`) |> 
  tidyr::gather(key = "year", value = "mda_binary", -subcounty) |> 
  mutate(year = as.integer(year),
         subcounty = as.character(subcounty)) |> 
  arrange(subcounty, year) |>                    
  group_by(subcounty) |>                         
  mutate(cum_mda = cumsum(mda_binary)) |> 
  ungroup()

# Add 1 year to mda before joining to take into account that MDA
# always occurs after the survey
mda_cum$year <- mda_cum$year + 1

sth <- left_join(sth, mda_cum)
sth$cum_mda[sth$year <= 2013] <- 0

sth_clean <- sth |> 
  st_drop_geometry() |> 
  select(prevASC, prevHKW, prevTT, cum_mda) |> 
  tidyr::gather(key = "specie", value = "prev", -cum_mda) |> 
  mutate(specie = factor(specie, levels = c("prevASC", "prevHKW", "prevTT"),
                         labels = c("Ascaris", "Hookworw", "Trichuris")))

# Prevalence vs. cum mda rounds
ggplot(sth_clean, aes(x = factor(cum_mda), y = prev)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.6) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.3) +
  facet_wrap(~ specie) +
  scale_y_log10(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Cumulative rounds of MDA",
    y = "STH prevalence (log scale)"
  ) +
  theme_minimal(base_size = 12)


ggsave(filename = "figures/cum_mda_kenya.pdf",
       device = cairo_pdf, 
       bg = "white", width = 10, height = 5)

# FIGURE 9. --------------------------------------------------------------------

# FIGURE 10. -------------------------------------------------------------------

# FIGURE 11. -------------------------------------------------------------------

## 1. IMPORT DATA --------------------------------------------------------------

# Survey data LF
lf <- readr::read_csv("data/lf_madagascar.csv")

# Madagascar boudaries
madagascar_bg <- st_read("data/mdg_adm0.shp")

# Madagascar IU boundaries
mdg_ius <- st_read("data/mdg_adm2.shp")

## 2. DATA PROCESSING ----------------------------------------------------------

# Set CRS
lf <- st_as_sf(lf, coords = c("lon", "lat"), crs = 4326)

# Define time panels LF
panel_labs_lf <- c("2004 (Baseline)", "2005–2008", "2009–2012",
                   "2013–2016", "2017–2020")

lf_pan <- lf |>
  mutate(
    panel = case_when(
      year == 2004 ~ "2004 (Baseline)",
      year %in% 2005:2008 ~ "2005–2008",
      year %in% 2009:2012 ~ "2009–2012",
      year %in% 2013:2016 ~ "2013–2016",
      year %in% 2017:2020 ~ "2017–2020",
      .default = NA_character_
    ),
    panel = factor(panel, levels = panel_labs_lf)
  )


base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "top"
  )

# A function to plot a single year-range panel
plot_panel <- function(lbl) {
  pdat <- filter(lf_pan, panel == lbl)
  ggplot() +
    geom_sf(data = madagascar_bg, fill = "grey95", color = "black") +
    geom_sf(
      data = pdat,
      aes(fill = prevalence),
      stroke = 0.2, shape = 21, size = 1.3, alpha = 1
    ) +
    scale_fill_distiller(
      type = "div",
      palette = "RdYlBu",
      direction = -1,
      limits = c(0, 0.75),
      breaks = seq(0, 0.75, by = 0.1),
      labels = scales::label_percent(accuracy = 1),
      name = "Prevalence (%)",
      guide = guide_colorbar(
        title.position = "top",        
        title.hjust = 0.5,             
        barwidth = unit(6, "in"),      
        barheight = unit(0.4, "cm")
      )) +
    coord_sf() +
    scale_x_continuous(breaks = seq(43, 51, by = 2)) +
    labs(title = lbl) +
    base_theme
}

# Build each panel
p_2004 <- plot_panel(panel_labs_lf[1])
p_0508 <- plot_panel(panel_labs_lf[2])
p_0912 <- plot_panel(panel_labs_lf[3])
p_1316 <- plot_panel(panel_labs_lf[4])
p_1720 <- plot_panel(panel_labs_lf[5])

# Subcounty-only panel
p_subc <- ggplot() +
  geom_sf(data = madagascar_bg, fill = "grey95", color = NA) +
  geom_sf(data = mdg_ius, fill = "grey", color = "grey40", linewidth = 0.15) +
  scale_x_continuous(breaks = seq(43, 51, by = 2)) + 
  labs(title = "Implementation units") +
  base_theme +
  theme(legend.position = "none") +
  annotation_scale(location = "br", width_hint = 0.25) +  
  annotation_north_arrow(location = "tl", which_north = "true", 
                         height = unit(.5, "cm"), width = unit(.5, "cm")) 

# Assemble 3×2 layout and add title/caption
plt <- (guide_area() / (
  p_2004 | p_0508 | p_0912
) /
  (p_1316 | p_1720 | p_subc)) +
  plot_layout(guides = "collect", heights = unit(c(1, 1, 1), c("cm", "null", "null"))) &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title = element_text(size = 11, face = "bold", vjust = 0.3),
    legend.text = element_text(size = 9),
  )

plt
ggsave(filename = "figures/lf_madagascar_panels.pdf",
       plot = plt, device = cairo_pdf, 
       bg = "white", width = 6.5, height = 7.5)

# FIGURE 12. -------------------------------------------------------------------


# FIGURE 13. -------------------------------------------------------------------

ggplot(lf, aes(x = factor(mda_rounds), y = prevalence)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  scale_y_log10(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Cumulative rounds of MDA",
    y = "LF prevalence (log scale)"
  ) +
  theme_minimal(base_size = 12)

ggsave(filename = "figures/cum_mda_mdg.pdf",
       device = cairo_pdf, 
       bg = "white", width = 8, height = 6)

# FIGURE 14. -------------------------------------------------------------------



