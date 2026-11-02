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
  pts_sf %>% mutate(lek_id = row_i$lek_id, date = row_i$date)
})

names(pts_by_lek_year) <- paste(files_tbl$lek_id, files_tbl$date, sep = "__")

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
  }
}
