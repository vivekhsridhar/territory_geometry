rm(list = ls())

## Set working directory to location of this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Load libraries and custom functions
library(sf)
library(purrr)
library(dplyr)
library(spatstat.geom)
library(spatstat.explore)
library(tidyverse)
library(zoo)

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

## Lek configuration table
lek_configs <- tibble(lek_id   = c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC"),
                      location = c("Velavadar", "Velavadar", "TalChhapar"),
                      suffix   = c("LEK1", "LEK2", "TC"),
                      shp_file = c("Velavadar_Lek1_Area.shp", "Velavadar_Lek2_Area.shp", "TalChhapar_Area.shp"))

## Output folder
out_dir <- file.path(script_dir, "processed_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Build master table of all files across all leks
files_tbl <- map_dfr(seq_len(nrow(lek_configs)), function(i) {
  
  cfg <- lek_configs[i, ]
  
  data_dirs <- list.dirs(file.path(root_dir, cfg$location), recursive = FALSE, full.names = TRUE)
  data_dirs <- data_dirs[grepl("_COORDINATES$", basename(data_dirs))]
  
  map_dfr(data_dirs, function(d) {
    
    data_label <- sub("_COORDINATES$", "", basename(d))
    csv_path <- list.files(d, pattern = paste0("_", cfg$suffix, "\\.csv$"), full.names = TRUE)
    
    if (length(csv_path) == 0) return(NULL)
    
    tibble(lek_id = cfg$lek_id, location = cfg$location, suffix = cfg$suffix, 
           shp_file = cfg$shp_file, data_label = data_label,
           date = parse_label_to_date(data_label, month_lookup), csv_path = csv_path[1])
  })
}) %>% arrange(lek_id, date)

## PCF controls (each lek x date)
n_r <- 200
s_max <- 4
pcf_hs <- 0.3
min_n <- 30
sigma_s <- 2.5
weight_cap_q <- 0.95
correction <- "translate"
min_neff <- 40

## MAIN
## Compute PCFs for all leks and all dates
curve_list <- list()
summary_list <- list()

for (i in seq_len(nrow(files_tbl))) {
  
  row <- files_tbl[i, ]
  message("Processing ", row$lek_id, " : ", row$data_label, " (", i, "/", nrow(files_tbl), ")")
  
  lek_polygon <- st_read(file.path(root_dir, row$location, row$shp_file), quiet = TRUE) |>
    st_transform(32643) |> st_zm(drop = TRUE)
  
  df <- read.csv(row$csv_path)
  n_pts <- nrow(df)
  lek_points <- st_as_sf(df, coords = c("pos_x", "pos_y"), crs = 32643)
  
  res <- compute_pcf(lek_polygon = lek_polygon, lek_points_sf = lek_points,
                     n_r = n_r, pcf_hs = pcf_hs, sigma_s = sigma_s, correction = correction)
  
  curve_list[[i]] <- res$g_df %>%
    mutate(lek_id = row$lek_id, data_label = row$data_label, date = row$date,
           n_points = n_pts, nn_median = res$nn_median, row$N_eff)
  
  summary_list[[i]] <- tibble(lek_id = row$lek_id, data_label = row$data_label,
                              date = row$date, n_points = n_pts, nn_median = res$nn_median, 
                              bw_sigma = as.numeric(res$sigma), bw_sigma_s = as.numeric(res$sigma_s), 
                              pcf_hs = res$pcf_hs, pcf_h = res$pcf_h)
}

pcf_curves <- bind_rows(curve_list)
pcf_summary <- bind_rows(summary_list)
pcf_curves <- pcf_curves[pcf_curves$n_points >= min_n,]
pcf_summary <- pcf_summary[pcf_summary$n_points >= min_n,]

## Write outputs
curves_out_file <- file.path(out_dir, "pcf_curve_ALL.csv")
summary_out_file <- file.path(out_dir, "pcf_summary_ALL.csv")
write.csv(pcf_curves, curves_out_file, row.names = FALSE)
write.csv(pcf_summary, summary_out_file, row.names = FALSE)

message("Saved curves to: ", curves_out_file)
message("Saved summary to: ", summary_out_file)

## Perform peak detection on computed PCF curves
## Peak detection parameters
lower_nnd_bound_s <- 0.5
min_prominence <- 0.02
min_peak_sep_s <- 0.5

## Apply peak detection to all lek × date PCFs
peak_table <- pcf_curves %>%
  group_by(lek_id, date) %>%
  group_modify(~{
    df <- .x %>% arrange(r)
    med_nnd <- unique(df$nn_median)
    n_pts   <- unique(df$n_points)
    
    peaks <- detect_peaks(df$s, df$g, med_nnd)
    
    # Attach Neff at the detected peak location, then filter by support.
    # detect_peaks() itself is unchanged.
    peaks <- peaks %>%
      mutate(N_eff = df$N_eff[match(r_peak, df$r)]) %>%
      filter(is.finite(N_eff), N_eff >= min_neff) %>%
      select(-N_eff)

    if (nrow(peaks) == 0) {
      return(tibble(s_peak = NA_real_,
                    r_peak = NA_real_,
                    g_peak = NA_real_,
                    peak_prominence = NA_real_,
                    peak_curvature = NA_real_,
                    n_peaks = 0L,
                    n_points = n_pts))
    }

    peaks %>% mutate(n_peaks = n(), n_points = n_pts)
  }) %>% ungroup()

## Save peak table
peaks_out_file <- file.path(out_dir, "pcf_peak_table_ALL.csv")
write.csv(peak_table, peaks_out_file, row.names = FALSE)

message("Saved peak table to: ", peaks_out_file)
