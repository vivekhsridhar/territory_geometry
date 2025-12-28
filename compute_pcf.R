rm(list = ls())

## Load libraries
library(sf)
library(purrr)
library(dplyr)
library(spatstat.geom)
library(spatstat.explore)

## Lookup table to convert file names to dates
month_lookup <- c(JAN = 1, FEB = 2, MAR = 3, APR = 4,
                  MAY = 5, JUN = 6, JUL = 7, AUG = 8,
                  SEP = 9, OCT = 10, NOV = 11, DEC = 12)

parse_label_to_date <- function(data_label, month_lookup) {
  yy <- as.integer(substr(data_label, 1, 2))
  mm <- unname(month_lookup[substr(data_label, 3, 5)])
  as.Date(sprintf("20%02d-%02d-01", yy, mm))
}

## Compute g_inhom(r)
compute_pcf <- function(lek_polygon, lek_points_sf, r_vals, r_min = 1.2, correction = "translate") {
  # Define observation window
  W <- as.owin(st_geometry(lek_polygon))
  
  # Create point pattern object
  pts <- st_coordinates(lek_points_sf)
  X <- ppp(pts[, 1], pts[, 2], window = W)
  stopifnot(all(inside.owin(X$x, X$y, W)))
  
  # Nearest-neighbour distances (for summary)
  nn <- nndist(X)
  nn_median <- median(nn)
  
  # Intensity estimate (inhomogeneous)
  sigma <- bw.ppl(X)
  lambda_hat <- density.ppp(X, sigma = sigma, edge = TRUE, at = "pixels")
  
  # Inhomogeneous pair correlation on a common r-grid
  g <- pcfinhom(X, lambda = lambda_hat, r = r_vals, correction = correction)
  
  # Use translate correction values by default
  g_df <- tibble(r = g$r, g = g$trans) %>% filter(r > r_min)
  
  list(g_df = g_df, nn_median = nn_median, sigma = sigma)
}

## Root code and data directory
root_dir <- "/Users/vivekhsridhar/Library/Mobile Documents/com~apple~CloudDocs/Documents/Data/SatelliteImagery/GoogleEarth"
setwd(root_dir)

## Lek configuration table
lek_configs <- tibble(lek_id   = c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC"),
  location = c("Velavadar", "Velavadar", "TalChhapar"), suffix   = c("LEK1", "LEK2", "TC"),
  shp_file = c("Velavadar_Lek1_Area.shp", "Velavadar_Lek2_Area.shp", "TalChhapar_Area.shp"))

## Comparability controls
scale_mode <- "global"   # global reference across all leks
r_max_mult <- 4          # r_max = r_max_mult * reference_median_nn
n_r <- 240               # number of r points in the curve (common across files)
r_min <- 1.2             # drop tiny-r behaviour (body length of a blackbuck)
correction <- "translate"

## Output folder
dir.create("derived_metrics", showWarnings = FALSE)

## Build master table of all files across all leks
files_tbl <- map_dfr(seq_len(nrow(lek_configs)), function(i) {
  cfg <- lek_configs[i, ]
  
  data_dirs <- list.dirs(file.path(root_dir, cfg$location), recursive = FALSE, full.names = TRUE)
  data_dirs <- data_dirs[grepl("_COORDINATES$", basename(data_dirs))]
  
  map_dfr(data_dirs, function(d) {
    data_label <- sub("_COORDINATES$", "", basename(d))
    csv_path <- list.files(d, pattern = paste0("_", cfg$suffix, "\\.csv$"), full.names = TRUE)
    
    if (length(csv_path) == 0) return(NULL)
    
    tibble(lek_id = cfg$lek_id, location = cfg$location, suffix = cfg$suffix, shp_file = cfg$shp_file, 
           data_label = data_label, date = parse_label_to_date(data_label, month_lookup), 
           csv_path = csv_path[1])
  })
}) %>% arrange(lek_id, date)

if (nrow(files_tbl) == 0) {
  stop("No CSV files found across leks.")
}

## First pass: compute a GLOBAL reference median NND for comparability
nn_pass <- map_dfr(seq_len(nrow(files_tbl)), function(i) {
  
  row <- files_tbl[i, ]
  
  lek_polygon <- st_read(file.path(root_dir, row$location, row$shp_file), quiet = TRUE) |> st_transform(32643) |> st_zm(drop = TRUE)
  
  df <- read.csv(row$csv_path)
  pts_sf <- st_as_sf(df, coords = c("pos_x", "pos_y"), crs = 32643)
  
  W <- as.owin(st_geometry(lek_polygon))
  xy <- st_coordinates(pts_sf)
  X <- ppp(xy[, 1], xy[, 2], window = W)
  
  tibble(nn_median = median(nndist(X)))
})

files_tbl <- bind_cols(files_tbl, nn_pass)

ref_median_nn <- median(files_tbl$nn_median, na.rm = TRUE)

r_max  <- r_max_mult * ref_median_nn
r_vals <- seq(0, r_max, length.out = n_r)

message("GLOBAL reference median NND = ", round(ref_median_nn, 2), " m")
message("Using GLOBAL r_max = ", round(r_max, 2), " m")
message("Using ", n_r, " r-values")

## Second pass: compute PCFs for all leks on the SAME r-grid
curve_list <- list()
summary_list <- list()

for (i in seq_len(nrow(files_tbl))) {
  
  row <- files_tbl[i, ]
  message("Processing ", row$lek_id, " : ", row$data_label, " (", i, "/", nrow(files_tbl), ")")
  
  lek_polygon <- st_read(file.path(root_dir, row$location, row$shp_file), quiet = TRUE) |> st_transform(32643) |> st_zm(drop = TRUE)
  
  df <- read.csv(row$csv_path)
  lek_points <- st_as_sf(df, coords = c("pos_x", "pos_y"), crs = 32643)
  
  res <- compute_pcf(lek_polygon = lek_polygon, lek_points_sf = lek_points, 
                     r_vals = r_vals, r_min = r_min, correction = correction)
  
  res <- compute_pcf(lek_polygon = lek_polygon, lek_points_sf = lek_points, 
                     r_vals = r_vals, correction = correction)
  
  # Save full curve
  curve_list[[i]] <- res$g_df %>% mutate(lek_id = row$lek_id, data_label = row$data_label,
                                         date = row$date, r_max_used = r_max, ref_median_nn = ref_median_nn)
  
  # Save per-date summary
  summary_list[[i]] <- tibble(lek_id = row$lek_id, data_label = row$data_label, date = row$date, 
                              n_points = nrow(df), nn_median = res$nn_median, bw_sigma = as.numeric(res$sigma),
                              ref_median_nn = ref_median_nn, r_max_used = r_max)
}

pcf_curves  <- bind_rows(curve_list)
pcf_summary <- bind_rows(summary_list)

## Write outputs
write.csv(pcf_curves, file.path("derived_metrics", "pcf_curve_ALL.csv"), row.names = FALSE)
write.csv(pcf_summary, file.path("derived_metrics", "pcf_summary_ALL.csv"), row.names = FALSE) 

message("Saved curves to:  derived_metrics/pcf_curve_ALL.csv")
message("Saved summary to: derived_metrics/pcf_summary_ALL.csv")
