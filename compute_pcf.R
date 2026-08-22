rm(list = ls())

## Load libraries
library(sf)
library(purrr)
library(dplyr)
library(spatstat.geom)
library(spatstat.explore)
library(tidyverse)
library(zoo)

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
lek_configs <- tibble(
  lek_id   = c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC"),
  location = c("Velavadar", "Velavadar", "TalChhapar"),
  suffix   = c("LEK1", "LEK2", "TC"),
  shp_file = c("Velavadar_Lek1_Area.shp", "Velavadar_Lek2_Area.shp", "TalChhapar_Area.shp")
)

## PCF controls (each lek x date)
n_r <- 200
s_max <- 4
pcf_hs <- 0.3
min_n <- 20
sigma_s <- 2
weight_cap_q <- 0.95
correction <- "translate"
min_neff <- 35

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
    
    tibble(lek_id = cfg$lek_id,
           location = cfg$location,
           suffix = cfg$suffix, 
           shp_file = cfg$shp_file,
           data_label = data_label,
           date = parse_label_to_date(data_label, month_lookup),
           csv_path = csv_path[1])
  })
}) %>% arrange(lek_id, date)

## HELPER FUNCTIONS
## Effective number of unordered point pairs contributing to each PCF value
compute_neff <- function(X, lambda_used, r_vals, pcf_h) {

  # Pairwise distances, translation edge weights, and inverse-intensity weights
  D <- pairdist(X)
  E <- edge.Trans(X)
  W_lambda <- outer(1 / lambda_used, 1 / lambda_used, "*")

  # Keep each unordered pair once
  upper <- upper.tri(D)

  purrr::map_dbl(r_vals, function(target_r) {

    # Epanechnikov kernel used by pcfinhom()
    u <- (target_r - D) / pcf_h
    K <- ifelse(abs(u) <= 1, (3 / (4 * pcf_h)) * (1 - u^2), 0)

    # Pair contribution, up to constants common to all pairs at target_r
    contribution <- K * E * W_lambda
    w <- contribution[upper & is.finite(contribution) & contribution > 0]

    p <- w / sum(w)
    1 / sum(p^2)
  })
}

## Compute g_inhom(r) per lek x date
compute_pcf <- function(lek_polygon, lek_points_sf, n_r = 200, s_max = 4,
                        pcf_hs = 0.2, sigma_s = 2, weight_cap_q = 0.95, 
                        correction = "translate") {
  
  # Define observation window
  W <- as.owin(st_geometry(lek_polygon))
  
  # Create point pattern object
  pts <- st_coordinates(lek_points_sf)
  X <- ppp(pts[, 1], pts[, 2], window = W)
  stopifnot(all(inside.owin(X$x, X$y, W)))
  
  # Median nearest-neighbour distance for this point pattern
  nn <- nndist(X)
  nn_median <- median(nn)
  
  # Define r-range from the curve's own median NND
  s_vals <- seq(0, s_max, length.out = n_r)
  r_vals <- s_vals * nn_median
  
  # Intensity estimate (inhomogeneous)
  sigma <- sigma_s * nn_median
  lambda_hat <- density.ppp(X, sigma = sigma, edge = TRUE, at = "points", leaveoneout = TRUE)
  
  # Cap lambda weighting (this prevents extreme effects of outliers on the pcf curves)
  inv_lambda <- 1 / lambda_hat
  weight_cap <- quantile(inv_lambda, weight_cap_q)
  inv_lambda_used <- pmin(inv_lambda, weight_cap)
  lambda_used <- 1 / inv_lambda_used
  
  # Inhomogeneous pair correlation
  pcf_h <- pcf_hs * nn_median
  g <- pcfinhom(X, lambda = lambda_used, r = r_vals, h = pcf_h, correction = correction)
  
  # Extract translation- and weight-corrected PCF
  g_df <- tibble(r = g$r, s = g$r / nn_median, g = g$trans)

  # Effective pair support on exactly the same r-grid as the PCF
  N_eff <- compute_neff(X = X, lambda_used = lambda_used, r_vals = g$r, pcf_h = pcf_h)

  g_df <- g_df %>% mutate(N_eff = N_eff)
  
  weight_df <- tibble(x = X$x, y = X$y, lambda_hat = lambda_hat,
                      inv_lambda = inv_lambda, inv_lambda_used = inv_lambda_used, 
                      capped = inv_lambda > weight_cap)
  
  list(g_df = g_df, weight_df = weight_df, nn_median = nn_median, 
       sigma = sigma, sigma_s = sigma / nn_median,
       pcf_h = pcf_h, pcf_hs = pcf_hs, weight_cap = weight_cap,
       n_capped = sum(inv_lambda > weight_cap))
}

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

## PEAK DETECTION
## Peak detection parameters
lower_nnd_mult <- 0.8
smooth_k <- 5
min_prominence <- 0.02
min_sep_mult <- 0.1

## Peak detection function
detect_peaks <- function(s, g, med_nnd) {
  
  # Rescale distance to go back to metres
  r <- s * med_nnd
  
  # Apply lower cutoff in rescaled units
  keep <- s >= lower_nnd_mult
  s <- s[keep]
  r <- r[keep]
  g <- g[keep]
  
  # Smooth g for robust peak detection and estimation of peak properties
  g_s <- zoo::rollmean(g, k = smooth_k, fill = NA, align = "center")
  
  # Local maxima on the smoothed curve
  is_peak <- g_s > dplyr::lag(g_s) & g_s > dplyr::lead(g_s)
  peak_idx <- which(is_peak)
  
  # Candidate peaks. Peak height is taken from the smoothed curve.
  peaks <- tibble(idx = peak_idx,
                  s_peak = s[peak_idx],
                  r_peak = r[peak_idx],
                  g_peak = g_s[peak_idx])
  
  # Prominence window width in rescaled units
  win_half_width <- 0.75
  
  # Step size on rescaled axis for second-derivative curvature
  ds <- median(diff(s), na.rm = TRUE)
  
  # Estimate prominence and curvature from the smoothed curve
  peaks <- purrr::pmap_dfr(peaks, function(idx, s_peak, r_peak, g_peak) {
    
    left_limit_s  <- s_peak - win_half_width
    right_limit_s <- s_peak + win_half_width
    
    left_idx  <- which(s >= left_limit_s & s < s_peak)
    right_idx <- which(s > s_peak & s <= right_limit_s)
    
    # A peak requires support on both sides to define prominence
    left_min  <- min(g_s[left_idx], na.rm = TRUE)
    right_min <- min(g_s[right_idx], na.rm = TRUE)
    baseline  <- max(left_min, right_min)
    peak_prominence <- g_peak - baseline
    
    # Curvature: negative second derivative of smoothed g(s)
    gpp <- (g_s[idx + 1] - 2 * g_s[idx] + g_s[idx - 1]) / (ds^2)
    peak_curvature <- -gpp
    
    tibble(idx = idx,
           s_peak = s_peak,
           r_peak = r_peak,
           g_peak = g_peak,
           peak_prominence = peak_prominence,
           peak_curvature = peak_curvature)
  }) %>%
    filter(!is.na(peak_prominence), peak_prominence >= min_prominence)
  
  if (nrow(peaks) == 0) return(NULL)
  
  # Iteratively enforce minimum peak separation.
  # At each step, retain the most prominent remaining peak and remove
  # all other candidates within min_sep_mult on the rescaled distance axis.
  remaining <- peaks %>% arrange(desc(peak_prominence), desc(g_peak))
  selected <- list()
  
  while (nrow(remaining) > 0) {
    chosen <- remaining[1, ]
    selected[[length(selected) + 1]] <- chosen
    
    remaining <- remaining %>% filter(abs(s_peak - chosen$s_peak) >= min_sep_mult)
  }
  
  bind_rows(selected) %>% arrange(s_peak) %>% select(-idx)
}

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
