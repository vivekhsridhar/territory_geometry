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
n_sims <- 10
nnd_min_quantile <- 0.10  # Minimum spacing for constrained null: 10th percentile of observed current NND
proposal_batch_size <- 10000
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

  pts_sf %>% mutate(lek_id = row_i$lek_id, date = row_i$date,
                    point_id = seq_len(nrow(pts_sf)) - 1L)
})

names(pts_by_lek_year) <- paste(files_tbl$lek_id, files_tbl$date, sep = "__")

## Initialise output containers
observed_list <- list()
null_unconstrained_list <- list()
null_constrained_list <- list()

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

    # Keep current territories inside the previous observation's KDE core
    pts_curr_in_core <- subset_points_to_kde_core(pts_curr, kde_prev)

    # Classify current locations using their distance to the nearest previous territory
    obs_nnd <- nnd(pts_curr_in_core, pts_prev)
    reused <- obs_nnd <= reuse_thresh

    # Remove newly observed core territories while preserving surrounding neighbours
    new_point_ids <- pts_curr_in_core$point_id[!reused]
    pts_without_new <- pts_curr %>% filter(!point_id %in% new_point_ids)
    n_new <- length(new_point_ids)

    if (any(reused)) {
      # Calculate six-neighbour CV before and after the observed additions
      focal_ids <- pts_curr_in_core$point_id[reused]

      before <- compute_local_regularity(lek_polygon, pts_without_new) %>%
        mutate(point_id = pts_without_new$point_id) %>%
        select(point_id, six_neighbour_cv)

      after <- compute_local_regularity(lek_polygon, pts_curr) %>%
        mutate(point_id = pts_curr$point_id)

      observed <- after %>%
        filter(point_id %in% focal_ids) %>%
        select(point_id, six_neighbour_cv) %>%
        inner_join(before, by = 'point_id', suffix = c('_after', '_before')) %>%
        mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr, n_new = n_new)

      observed_list[[length(observed_list) + 1]] <- observed

      if (n_new > 0) {
        # Spacing threshold for the NND-constrained null
        xy_curr <- st_coordinates(pts_curr)
        X_curr <- ppp(xy_curr[, 1], xy_curr[, 2],
                      window = as.owin(st_geometry(lek_polygon)))
        min_spacing <- unname(quantile(nndist(X_curr), nnd_min_quantile))

        # Both nulls sample from the previous KDE, restricted to its core
        proposal_intensity <- kde_prev$intensity
        intensity_values <- proposal_intensity$v
        cutoff <- kde_prev$contour_level * max(intensity_values, na.rm = TRUE)
        proposal_intensity$v[!is.finite(intensity_values) | intensity_values < cutoff] <- 0

        # Null 1: KDE-weighted placement without a current-year spacing constraint
        pool <- sample_new_locations(pts_prev, lek_polygon, proposal_intensity,
                                     n_new * n_sims, reuse_thresh,
                                     proposal_batch_size)

        for (sim in seq_len(n_sims) - 1L) {
          # Draw exactly n_new locations under each null
          start <- sim * n_new + 1L
          end <- start + n_new - 1L
          
          pts_new <- pool$points[start:end, ]
          pts_new$point_id <- new_point_ids
          pts_sim <- rbind(pts_without_new %>% select(point_id),
                           pts_new %>% select(point_id))

          null_unconstrained_list[[length(null_unconstrained_list) + 1]] <-
            compute_local_regularity(lek_polygon, pts_sim) %>%
            mutate(point_id = pts_sim$point_id) %>%
            filter(point_id %in% focal_ids) %>%
            select(point_id, six_neighbour_cv) %>%
            mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr,
                   sim = sim, n_new = n_new)

          # Null 2: new locations separated from fixed current points and each other
          draw <- sample_new_locations_spaced(
            pts_prev, pts_without_new, lek_polygon, proposal_intensity,
            n_new, reuse_thresh, min_spacing
          )

          pts_new <- draw$points
          pts_new$point_id <- new_point_ids
          pts_sim <- rbind(pts_without_new %>% select(point_id),
                           pts_new %>% select(point_id))

          null_constrained_list[[length(null_constrained_list) + 1]] <-
            compute_local_regularity(lek_polygon, pts_sim) %>%
            mutate(point_id = pts_sim$point_id) %>%
            filter(point_id %in% focal_ids) %>%
            select(point_id, six_neighbour_cv) %>%
            mutate(lek_id = lek, date_prev = date_prev, date_curr = date_curr,
                   sim = sim, n_new = n_new, min_spacing = min_spacing)
        }
      }
    }
  }
}

## Export territory-level measurements only; comparisons belong in the notebook
out_dir <- file.path(getwd(), 'processed_data')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

observed_tbl <- bind_rows(observed_list) %>% arrange(lek_id, date_curr, point_id)
unconstrained_tbl <- bind_rows(null_unconstrained_list) %>% arrange(lek_id, date_curr, sim, point_id)
constrained_tbl <- bind_rows(null_constrained_list) %>% arrange(lek_id, date_curr, sim, point_id)

write_csv(observed_tbl, file.path(out_dir, 'spatial_addition_observed.csv'))
write_csv(unconstrained_tbl, file.path(out_dir, 'spatial_addition_null_unconstrained.csv'))
write_csv(constrained_tbl, file.path(out_dir, 'spatial_addition_null_nnd_constrained.csv'))
