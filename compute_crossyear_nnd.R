rm(list = ls())

## Set working directory to location of this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Load libraries and custom functions
library(sf)
library(rstudioapi)
library(purrr)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(spatstat.geom)
library(spatstat.explore)

source('spatial_analysis_functions.R')

## SETUP / HOUSEKEEPING
## Get directory where this script lives
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

## Lookup table to convert file names to dates
month_lookup <- c(JAN = 1, FEB = 2, MAR = 3, APR = 4,
                  MAY = 5, JUN = 6, JUL = 7, AUG = 8,
                  SEP = 9, OCT = 10, NOV = 11, DEC = 12)

parse_label_to_date <- function(data_label, month_lookup) {
  yy <- as.integer(substr(data_label, 1, 2))
  mm <- unname(month_lookup[substr(data_label, 3, 5)])
  as.Date(sprintf("20%02d-%02d-01", yy, mm))
}

## Root code and data directory
root_dir <- "/Users/vivekhsridhar/Library/Mobile Documents/com~apple~CloudDocs/Documents/Data/SatelliteImagery/GoogleEarth"
setwd(root_dir)

## Output folders
out_dir <- file.path(script_dir, "processed_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir <- file.path(script_dir, "output")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

## Lek configurations
lek_configs <- tibble(lek_id   = c("Tal Chhapar", "Velavadar Lek 1", "Velavadar Lek 2"),
                      location = c("TalChhapar", "Velavadar", "Velavadar"),
                      suffix   = c("TC", "LEK1", "LEK2"),
                      shp_file = c("TalChhapar_Area.shp", "Velavadar_Lek1_Area.shp", "Velavadar_Lek2_Area.shp"))

plot_crossyear = FALSE

## Build master table of all files across all leks
files_tbl <- map_dfr(seq_len(nrow(lek_configs)), function(i) {
  
  cfg <- lek_configs %>% slice(i)
  
  data_dirs <- list.dirs(file.path(root_dir, cfg$location), recursive = FALSE, full.names = TRUE)
  data_dirs <- data_dirs[grepl("_COORDINATES$", basename(data_dirs))]
  
  map_dfr(data_dirs, function(d) {
    
    data_label <- sub("_COORDINATES$", "", basename(d))
    csv_path <- list.files(d, pattern = paste0("_", cfg$suffix, "\\.csv$"), full.names = TRUE)
    
    if (length(csv_path) == 0) return(NULL)
    
    tibble(lek_id = cfg$lek_id,
           location = cfg$location,
           suffix = cfg$suffix,
           shp_file = cfg$shp_file,
           data_label = data_label,
           date = parse_label_to_date(data_label, month_lookup),
           csv_path = csv_path[1])
  })
}) %>% arrange(lek_id, date)

## ANALYSIS PARAMETERS
## KDE settings
core_prob <- 0.8
kde_dimyx <- 512
min_points_kde <- 10

## Null simulation settings
n_sims <- 1000
set.seed(123)

## MAIN
## Read and store information about the polygons for each lek
lek_polygons <- map(seq_len(nrow(lek_configs)), function(i) {
  cfg <- lek_configs %>% slice(i)
  poly_sf <- st_read(file.path(root_dir, cfg$location, cfg$shp_file),
                     quiet = TRUE) %>% st_transform(32643) %>% st_zm(drop = TRUE) %>% st_make_valid()
})

names(lek_polygons) <- lek_configs$lek_id

## Compute a single representative KDE bandwidth per lek
sigma_fixed_tbl <- map_dfr(seq_len(nrow(lek_configs)), function(i) {
  # Extract config for this lek and get all files ordered by date
  cfg <- lek_configs %>% slice(i)
  files_lek <- files_tbl %>% filter(lek_id == cfg$lek_id) %>% arrange(date)
  
  lek_polygon <- lek_polygons[[cfg$lek_id]]
  
  # Loop over dates and estimate the bandwidth value for each time point
  sigma_tbl <- map_dfr(seq_len(nrow(files_lek)), function(j) {
    
    row_j <- files_lek[j, ]
    
    # Read the point coordinates and only keep ones inside the lek polygon
    pts_sf <- read.csv(row_j$csv_path) %>% st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
    pts_sf <- clip_points_to_polygon(pts_sf, lek_polygon)
    
    # Store the number of points and the KDE bandwidth 
    tibble(date = row_j$date, n_points = nrow(pts_sf),
           sigma_year = if (nrow(pts_sf) >= min_points_kde) get_kde_sigma(pts_sf, lek_polygon) else NA_real_)
  })
  
  # Take the median bandwidth across years as the fixed value for each lek
  tibble(lek_id = cfg$lek_id, sigma_fixed = median(sigma_tbl$sigma_year, na.rm = TRUE))
})

## Read and store all points (lek id x year)
pts_by_lek_year <- map(seq_len(nrow(files_tbl)), function(i) {
  
  row_i <- files_tbl[i, ]
  lek_polygon <- lek_polygons[[row_i$lek_id]]
  
  pts_sf <- read.csv(row_i$csv_path) %>% st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
  pts_sf <- clip_points_to_polygon(pts_sf, lek_polygon)
  pts_sf %>% mutate(lek_id = row_i$lek_id, date = row_i$date)
})

names(pts_by_lek_year) <- paste(files_tbl$lek_id, files_tbl$date, sep = "__")

## Initialise output containers
summary_list <- list()
pointwise_list <- list()
sim_list <- list()
pointwise_randomisation_list <- list()

## Split the file table by lek
files_by_lek <- split(files_tbl, files_tbl$lek_id)

## Compute fixed plot limits around all points for each lek
plot_limits_by_lek <- map(names(files_by_lek), function(lek) {

  files_lek <- files_by_lek[[lek]]

  all_xy <- map_dfr(seq_len(nrow(files_lek)), function(i) {
    pts <- pts_by_lek_year[[paste(files_lek$lek_id[i], files_lek$date[i], sep = "__")]]
    xy <- st_coordinates(pts)
    tibble(x = xy[, 1], y = xy[, 2])
  })

  pad <- 20

  list(xlim = c(min(all_xy$x) - pad, max(all_xy$x) + pad),
       ylim = c(min(all_xy$y) - pad, max(all_xy$y) + pad))
})

names(plot_limits_by_lek) <- names(files_by_lek)

## For each lek
for (lek in names(files_by_lek)) {
  
  # Get all files ordered by date
  files_lek <- files_by_lek[[lek]] %>% arrange(date)
  
  # Retreive the lek polygon and the KDE bandwidth for this lek
  lek_polygon <- lek_polygons[[lek]]
  sigma_fixed <- as.numeric(sigma_fixed_tbl[sigma_fixed_tbl$lek_id == lek,'sigma_fixed'])
  
  # Loop over consecutive dates
  for (i in 2:nrow(files_lek)) {
    
    # Extract points and metadata for previous and current year
    row_prev <- files_lek[i-1,]
    row_curr <- files_lek[i,]
    
    date_prev <- row_prev$date
    date_curr <- row_curr$date
    
    pts_prev <- pts_by_lek_year[[paste(row_prev$lek_id, row_prev$date, sep = "__")]]
    pts_curr <- pts_by_lek_year[[paste(row_curr$lek_id, row_curr$date, sep = "__")]]
    
    # Compute KDE grid from t-1
    kde_prev <- make_kde_grid(pts_sf = pts_prev, lek_polygon = lek_polygon, sigma = sigma_fixed, core_prob = core_prob, dimyx = kde_dimyx)
    
    # Keep only current-year points inside previous-year KDE core
    pts_curr_in_core <- subset_points_to_kde_core(pts_curr, kde_prev)
    
    # Plot current-year point pattern against the previous-year KDE core
    if (plot_crossyear == TRUE) {
      plot_crossyear_core_overlap(kde_grid = kde_prev, lek_polygon = lek_polygon,
                                  pts_prev = pts_prev, pts_curr = pts_curr,
                                  pts_curr_in_core = pts_curr_in_core,
                                  lek = lek, date_prev = date_prev,
                                  date_curr = date_curr, plot_dir = plot_dir,
                                  plot_limits = plot_limits_by_lek[[lek]])
    }
    
    # Compute nearest-neighbour distance from current to previous points
    obs_nnd <- nnd(pts_curr_in_core, pts_prev)
    dist_from_centre <- distance_to_kde_mode(pts_curr_in_core, kde_prev, crs_use = st_crs(pts_prev))
    
    # Simulate transformed previous-year points and compute NNDs
    sim_out <- simulate_transform_crossyear_nnd(prev_pts = pts_prev, curr_pts = pts_curr_in_core, n_sims = n_sims)
    sim_tbl <- sim_out$summary
    sim_pointwise_tbl <- sim_out$pointwise
    
    # Observed summary statistics
    obs_mean <- mean(obs_nnd)
    obs_median <- median(obs_nnd)
    
    # Qunatiles of simulated transformed mean NNDs
    rand_mean_transform_q025 <- quantile(sim_tbl$mean_nnd_transform_rand, 0.025, na.rm = TRUE)
    rand_mean_transform_q25  <- quantile(sim_tbl$mean_nnd_transform_rand, 0.25, na.rm = TRUE)
    rand_mean_transform_q50  <- quantile(sim_tbl$mean_nnd_transform_rand, 0.50, na.rm = TRUE)
    rand_mean_transform_q75  <- quantile(sim_tbl$mean_nnd_transform_rand, 0.75, na.rm = TRUE)
    rand_mean_transform_q975 <- quantile(sim_tbl$mean_nnd_transform_rand, 0.975, na.rm = TRUE)

    # Quantiles of simulated transformed median NNDs
    rand_median_transform_q025 <- quantile(sim_tbl$median_nnd_transform_rand, 0.025, na.rm = TRUE)
    rand_median_transform_q25  <- quantile(sim_tbl$median_nnd_transform_rand, 0.25, na.rm = TRUE)
    rand_median_transform_q50  <- quantile(sim_tbl$median_nnd_transform_rand, 0.50, na.rm = TRUE)
    rand_median_transform_q75  <- quantile(sim_tbl$median_nnd_transform_rand, 0.75, na.rm = TRUE)
    rand_median_transform_q975 <- quantile(sim_tbl$median_nnd_transform_rand, 0.975, na.rm = TRUE)

    # Store summary statistics for this time step
    summary_list[[length(summary_list) + 1]] <- tibble(lek_id = lek, date_prev = date_prev, date_curr = date_curr,
                                                       n_prev = nrow(pts_prev), n_curr = nrow(pts_curr), 
                                                       n_curr_in_prev_core = nrow(pts_curr_in_core),
                                                       obs_mean_nnd = obs_mean, obs_median_nnd = obs_median,
                                                       rand_mean_transform_nnd_q025 = as.numeric(rand_mean_transform_q025),
                                                       rand_mean_transform_nnd_q25 = as.numeric(rand_mean_transform_q25),
                                                       rand_mean_transform_nnd_q50 = as.numeric(rand_mean_transform_q50),
                                                       rand_mean_transform_nnd_q75 = as.numeric(rand_mean_transform_q75),
                                                       rand_mean_transform_nnd_q975 = as.numeric(rand_mean_transform_q975),
                                                       rand_median_transform_nnd_q025 = as.numeric(rand_median_transform_q025),
                                                       rand_median_transform_nnd_q25 = as.numeric(rand_median_transform_q25),
                                                       rand_median_transform_nnd_q50 = as.numeric(rand_median_transform_q50),
                                                       rand_median_transform_nnd_q75 = as.numeric(rand_median_transform_q75),
                                                       rand_median_transform_nnd_q975 = as.numeric(rand_median_transform_q975))
    
    # Store point-level NNDs for current points within the KDE core from previous year's points
    pointwise_list[[length(pointwise_list) + 1]] <- pts_curr_in_core %>%
      mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr, point_id = seq_along(obs_nnd),
             nnd_to_prev = obs_nnd, dist_from_centre = dist_from_centre) %>% st_drop_geometry()
    
    # Store simulation output for the same transition
    sim_list[[length(sim_list) + 1]] <- sim_tbl %>%
      mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr)

    # Store point-level simulation output for the same transition
    pointwise_randomisation_list[[length(pointwise_randomisation_list) + 1]] <- sim_pointwise_tbl %>%
      mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr,
             dist_from_centre = dist_from_centre[point_id])
  }
}

summary_tbl <- bind_rows(summary_list) %>% arrange(lek_id, date_curr)
pointwise_tbl <- bind_rows(pointwise_list) %>% arrange(lek_id, date_curr)
sim_tbl_all <- bind_rows(sim_list) %>% arrange(lek_id, date_curr, sim)
pointwise_randomisation_tbl <- bind_rows(pointwise_randomisation_list) %>% arrange(lek_id, date_curr, sim, point_id)

## Export dataframes
summary_tbl <- summary_tbl %>%
  transmute(lek_id, date_prev, date_now = date_curr, n_prev = n_prev, n_curr = n_curr, n_curr_in_prev_core = n_curr_in_prev_core,
            mean_crossyear_nnd = obs_mean_nnd, median_crossyear_nnd = obs_median_nnd)

simulation_tbl <- sim_tbl_all %>%
  transmute(lek_id, date_prev, date_now = date_curr, sim, 
            mean_crossyear_nnd_transform_rand = mean_nnd_transform_rand,
            median_crossyear_nnd_transform_rand = median_nnd_transform_rand)

pointwise_randomisation_tbl <- pointwise_randomisation_tbl %>%
  transmute(lek_id, date_prev, date_now = date_curr, sim, point_id, nnd_to_prev, dist_from_centre)

crossyear_nnd_tbl <- simulation_tbl %>%
  left_join(summary_tbl, by = c("lek_id", "date_prev", "date_now")) %>%
  relocate(lek_id, date_prev, date_now, n_curr_in_prev_core, sim, mean_crossyear_nnd, median_crossyear_nnd,
           mean_crossyear_nnd_transform_rand, median_crossyear_nnd_transform_rand)

write_csv(summary_tbl, file.path(out_dir, "crossyear_nnd_summary.csv"))
write_csv(pointwise_tbl, file.path(out_dir, "crossyear_nnd_pointwise.csv"))
write_csv(pointwise_randomisation_tbl, file.path(out_dir, "crossyear_nnd_pointwise_randomisations.csv"))
write_csv(crossyear_nnd_tbl, file.path(out_dir, "crossyear_nnd_with_randomisations.csv"))
