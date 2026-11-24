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

## MAIN
## Read lek polygons
lek_polygons <- map(seq_len(nrow(lek_configs)), function(i) {
  cfg <- lek_configs %>% slice(i)
  st_read(file.path(root_dir, cfg$location, cfg$shp_file), quiet = TRUE) %>%
    st_transform(32643) %>% st_zm(drop = TRUE) %>% st_make_valid()
})
names(lek_polygons) <- lek_configs$lek_id

## Read each observation and attach six-neighbour CV computed from all its territories
pts_by_lek_year <- map(seq_len(nrow(files_tbl)), function(i) {
  row_i <- files_tbl[i, ]
  lek_polygon <- lek_polygons[[row_i$lek_id]]
  pts_sf <- read.csv(row_i$csv_path) %>%
    st_as_sf(coords = c("pos_x", "pos_y"), crs = 32643)
  pts_sf <- clip_points_to_polygon(pts_sf, lek_polygon)
  local_metrics <- compute_local_regularity(lek_polygon, pts_sf)
  pts_sf %>%
    mutate(lek_id = row_i$lek_id, date = row_i$date,
           point_id = seq_len(nrow(pts_sf)) - 1L) %>%
    bind_cols(local_metrics)
})
names(pts_by_lek_year) <- paste(files_tbl$lek_id, files_tbl$date, sep = "__")

retention_list <- list()
files_by_lek <- split(files_tbl, files_tbl$lek_id)

for (lek in names(files_by_lek)) {
  files_lek <- files_by_lek[[lek]] %>% arrange(date)
  lek_polygon <- lek_polygons[[lek]]

  for (i in 2:nrow(files_lek)) {
    row_prev <- files_lek[i - 1, ]
    row_curr <- files_lek[i, ]
    date_prev <- row_prev$date
    date_curr <- row_curr$date

    pts_prev <- pts_by_lek_year[[paste(row_prev$lek_id, date_prev, sep = "__")]]
    pts_curr <- pts_by_lek_year[[paste(row_curr$lek_id, date_curr, sep = "__")]]

    # Select previous territories inside the previous KDE core
    sigma_prev <- get_kde_sigma(pts_prev, lek_polygon)
    kde_prev <- make_kde_grid(pts_prev, lek_polygon, sigma_prev,
                             core_thresh = core_thresh, dimyx = kde_dimyx)
    pts_prev_in_core <- subset_points_to_kde_core(pts_prev, kde_prev)

    # Did each previous core location persist into the current observation?
    nnd_to_next <- nnd(pts_prev_in_core, pts_curr)
    retention_list[[length(retention_list) + 1]] <- pts_prev_in_core %>%
      mutate(date_prev = date_prev, date_curr = date_curr,
             interval_days = as.numeric(date_curr - date_prev),
             nnd_to_next = nnd_to_next,
             persistent = nnd_to_next <= reuse_thresh) %>%
      st_drop_geometry()
  }
}

## Export territory-level results; preserve existing filename for downstream analyses
retention_tbl <- bind_rows(retention_list) %>% arrange(lek_id, date_curr, point_id)
out_dir <- file.path(getwd(), 'processed_data')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(retention_tbl, file.path(out_dir, 'selective_retention.csv'))
