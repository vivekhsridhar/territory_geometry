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

## Reuse threshold derived from the drone marking measurements
sizes <- read_tsv('rawdata/territory_sizes_2023.txt', show_col_types = FALSE)
reuse_thresh <- unname(quantile(sizes$sizes, 0.95))

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
root_dir <- file.path(getwd(), 'rawdata')

## Lek configurations
lek_configs <- tibble(lek_id   = c("Tal Chhapar", "Velavadar Lek 1", "Velavadar Lek 2"),
                      location = c("TalChhapar", "Velavadar", "Velavadar"),
                      suffix   = c("TC", "LEK1", "LEK2"),
                      shp_file = c("TalChhapar_Area.shp", "Velavadar_Lek1_Area.shp", "Velavadar_Lek2_Area.shp"))

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

## KDE settings
core_thresh <- 0.2
kde_dimyx <- 512

## Null simulation settings
n_sims <- 1000
n_preview_sims <- 3
proposal_batch_size <- 10000
max_proposals_per_point <- 1000
set.seed(123)

## MAIN
## Read and store information about the polygons for each lek
lek_polygons <- map(seq_len(nrow(lek_configs)), function(i) {
  cfg <- lek_configs %>% slice(i)
  poly_sf <- st_read(file.path(root_dir, cfg$location, cfg$shp_file),
                     quiet = TRUE) %>% st_transform(32643) %>% st_zm(drop = TRUE) %>% st_make_valid()
})

names(lek_polygons) <- lek_configs$lek_id

## Read and store all points (lek id x year)
pts_by_lek_year <- map(seq_len(nrow(files_tbl)), function(i) {
  
  row_i <- files_tbl[i, ]
  lek_polygon <- lek_polygons[[row_i$lek_id]]
  
  pts_sf <- read.csv(row_i$csv_path) %>% st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
  pts_sf <- clip_points_to_polygon(pts_sf, lek_polygon)

  # Attach local spacing measurements to each territory
  local_metrics <- compute_local_regularity(lek_polygon, pts_sf)
  pts_sf %>% mutate(lek_id = row_i$lek_id, date = row_i$date, point_id = seq_len(nrow(pts_sf)) - 1L) %>%
    bind_cols(local_metrics)
})

names(pts_by_lek_year) <- paste(files_tbl$lek_id, files_tbl$date, sep = "__")

## Initialise output containers
retention_list <- list()
pointwise_list <- list()
addition_list <- list()
addition_summary_list <- list()
sim_list <- list()
sim_pointwise_list <- list()
core_grid_list <- list()
summary_list <- list()

## Split the file table by lek
files_by_lek <- split(files_tbl, files_tbl$lek_id)

## For each lek
for (lek in names(files_by_lek)) {
  
  # Get all files ordered by date
  files_lek <- files_by_lek[[lek]] %>% arrange(date)
  
  # Retrieve the lek polygon
  lek_polygon <- lek_polygons[[lek]]
  
  # Loop over the next available observations
  for (i in 2:nrow(files_lek)) {
    
    # Extract points and metadata for the previous and current observations
    row_prev <- files_lek[i-1,]
    row_curr <- files_lek[i,]
    
    date_prev <- row_prev$date
    date_curr <- row_curr$date
    
    pts_prev <- pts_by_lek_year[[paste(row_prev$lek_id, row_prev$date, sep = "__")]]
    pts_curr <- pts_by_lek_year[[paste(row_curr$lek_id, row_curr$date, sep = "__")]]
    
    sigma_prev <- get_kde_sigma(pts_prev, lek_polygon)
      
    # Compute the KDE grid from the previous observation
    kde_prev <- make_kde_grid(pts_sf = pts_prev, lek_polygon = lek_polygon, sigma = sigma_prev, core_thresh = core_thresh, dimyx = kde_dimyx)

    # Store the core intensity grid for spatial diagnostics
    core_grid_list[[length(core_grid_list) + 1]] <- kde_prev$kde_df %>% filter(in_core) %>% select(x, y, p) %>%
      mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr, xstep = kde_prev$xstep, ystep = kde_prev$ystep)
    
    # Keep current territories inside the previous observation's KDE core
    pts_curr_in_core <- subset_points_to_kde_core(pts_curr, kde_prev)

    # Track persistence of previous territories inside their KDE core
    pts_prev_in_core <- subset_points_to_kde_core(pts_prev, kde_prev)
    nnd_to_next <- nnd(pts_prev_in_core, pts_curr)
    retention_list[[length(retention_list) + 1]] <- pts_prev_in_core %>%
      mutate(date_prev = date_prev, date_curr = date_curr, interval_days = as.numeric(date_curr - date_prev),
             nnd_to_next = nnd_to_next, persistent = nnd_to_next <= reuse_thresh,
             persistence_edge_ok = boundary_distance >= reuse_thresh) %>% st_drop_geometry()

    # Classify current locations using their distance to the nearest previous territory
    obs_nnd <- nnd(pts_curr_in_core, pts_prev)
    reused <- obs_nnd <= reuse_thresh

    # Store current territories and their spatial reuse classification
    pts_curr_in_core <- pts_curr_in_core %>% mutate(nnd_to_prev = obs_nnd, reused = reused)
    pointwise_list[[length(pointwise_list) + 1]] <- pts_curr_in_core %>%
      mutate(date_prev = date_prev, date_curr = date_curr) %>% st_drop_geometry()

    # Remove newly observed core territories while preserving surrounding neighbours
    nnd_reference <- median(pts_curr$nnd)
    new_point_ids <- pts_curr_in_core$point_id[!reused]
    pts_without_new <- pts_curr %>% filter(!point_id %in% new_point_ids)
    n_new <- length(new_point_ids)
    n_null_requested <- 0L
    n_null_complete <- 0L
    null_n_proposals <- 0L
    null_acceptance_rate <- NA_real_

    if (any(reused)) {
      # Measure spacing before additions using the fixed current-observation reference
      before <- compute_local_regularity(lek_polygon, pts_without_new, nnd_reference) %>%
        mutate(point_id = pts_without_new$point_id)

      # Compare the same reused focal territories before and after additions
      focal_ids <- pts_curr_in_core$point_id[reused]
      observed_comparison <- compare_addition_regularity(lek_polygon, pts_curr, before, focal_ids, nnd_reference)
      addition_list[[length(addition_list) + 1]] <- observed_comparison %>%
        mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr,
               n_reused = length(focal_ids), n_new = n_new)

      observed_summary <- summarise_addition_regularity(observed_comparison)
      addition_summary_list[[length(addition_summary_list) + 1]] <- observed_summary %>%
        mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr, n_new = n_new, nnd_reference = nnd_reference)

      if (n_new > 0 && any(observed_comparison$nnd_comparable | observed_comparison$six_neighbour_comparable)) {
        # Restrict the previous KDE intensity to the existing core
        proposal_intensity <- kde_prev$intensity
        intensity_values <- proposal_intensity$v
        cutoff <- kde_prev$contour_level * max(intensity_values, na.rm = TRUE)
        proposal_intensity$v[!is.finite(intensity_values) | intensity_values < cutoff] <- 0

        # Draw independent feasible locations for all null configurations
        n_null_requested <- n_sims
        pool <- sample_new_locations(pts_prev, lek_polygon, proposal_intensity, n_new * n_sims, reuse_thresh,
                                     proposal_batch_size, max_proposals_per_point)
        n_null_complete <- nrow(pool$points) %/% n_new
        null_n_proposals <- pool$n_proposals
        null_acceptance_rate <- pool$n_accepted / pool$n_proposals

        # Compare complete null configurations at the observed reused locations
        sim_results <- map(seq_len(n_sims) - 1L, function(s) {
          n_placed <- max(0, min(n_new, nrow(pool$points) - s * n_new))
          if (n_placed < n_new) return(list(
            summary = tibble(sim = s, complete = FALSE, n_placed = n_placed),
            pointwise = tibble(sim = integer(), point_id = integer(), pos_x = numeric(), pos_y = numeric())))

          pts_new <- pool$points[s * n_new + seq_len(n_new), ] %>%
            mutate(point_id = new_point_ids) %>% select(point_id)
          pts_sim <- rbind(pts_without_new %>% select(point_id), pts_new)
          comparison <- compare_addition_regularity(lek_polygon, pts_sim, before, focal_ids, nnd_reference)

          # Retain a few simulated configurations for plotting diagnostics
          preview <- NULL
          if (s < n_preview_sims) {
            xy <- st_coordinates(pts_new)
            preview <- tibble(sim = s, point_id = new_point_ids, pos_x = xy[, 1], pos_y = xy[, 2])
          }

          list(summary = summarise_addition_regularity(comparison) %>% mutate(sim = s, complete = TRUE, n_placed = n_placed),
               pointwise = preview)
        })

        sim_list[[length(sim_list) + 1]] <- map_dfr(sim_results, 'summary') %>%
          mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr, n_new = n_new, nnd_reference = nnd_reference)
        sim_pointwise_list[[length(sim_pointwise_list) + 1]] <- map_dfr(sim_results, 'pointwise') %>%
          mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr)
      }
    }

    # Record the counts defining each transition and its comparison
    summary_list[[length(summary_list) + 1]] <- tibble(lek_id = lek, date_prev = date_prev, date_curr = date_curr,
      interval_days = as.numeric(date_curr - date_prev), n_prev = nrow(pts_prev), n_prev_in_core = nrow(pts_prev_in_core),
      n_persistent = sum(nnd_to_next <= reuse_thresh), n_curr = nrow(pts_curr), n_curr_in_prev_core = nrow(pts_curr_in_core),
      n_reused = sum(reused), n_new = sum(!reused), n_background = nrow(pts_curr) - nrow(pts_curr_in_core),
      nnd_reference = nnd_reference, reuse_thresh = reuse_thresh,
      n_null_requested = n_null_requested, n_null_complete = n_null_complete,
      null_n_proposals = null_n_proposals, null_acceptance_rate = null_acceptance_rate)
  }
}

## Combine territory-level results across transitions
local_regularity_sf <- bind_rows(pts_by_lek_year)
local_xy <- st_coordinates(local_regularity_sf)
local_regularity_tbl <- local_regularity_sf %>%
  mutate(pos_x = local_xy[, 1], pos_y = local_xy[, 2]) %>% st_drop_geometry()
retention_tbl <- bind_rows(retention_list) %>% arrange(lek_id, date_curr, point_id)
pointwise_tbl <- bind_rows(pointwise_list) %>% arrange(lek_id, date_curr, point_id)
addition_tbl <- bind_rows(addition_list) %>% arrange(lek_id, date_curr, point_id)
addition_summary_tbl <- bind_rows(addition_summary_list) %>% arrange(lek_id, date_curr)
sim_tbl_all <- bind_rows(sim_list) %>% arrange(lek_id, date_curr, sim)
sim_pointwise_tbl <- bind_rows(sim_pointwise_list) %>% arrange(lek_id, date_curr, sim, point_id)
core_grid_tbl <- bind_rows(core_grid_list) %>% arrange(lek_id, date_curr)
summary_tbl <- bind_rows(summary_list) %>% arrange(lek_id, date_curr)

## Export dataframes
out_dir <- file.path(getwd(), 'processed_data')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(local_regularity_tbl, file.path(out_dir, 'reuse_local_regularity.csv'))
write_csv(retention_tbl, file.path(out_dir, 'reuse_retention_pointwise.csv'))
write_csv(pointwise_tbl, file.path(out_dir, 'reuse_current_pointwise.csv'))
write_csv(addition_tbl, file.path(out_dir, 'reuse_addition_pointwise.csv'))
write_csv(addition_summary_tbl, file.path(out_dir, 'reuse_addition_summary.csv'))
write_csv(sim_tbl_all, file.path(out_dir, 'reuse_addition_null.csv'))
write_csv(sim_pointwise_tbl, file.path(out_dir, 'reuse_addition_null_preview.csv'))
write_csv(core_grid_tbl, file.path(out_dir, 'reuse_kde_core.csv'))
write_csv(summary_tbl, file.path(out_dir, 'reuse_transition_summary.csv'))
