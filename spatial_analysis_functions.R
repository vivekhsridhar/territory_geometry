## HELPER FUNCTIONS
## 1) Compute NND metrics
compute_nnd <- function(lek_polygon, lek_points) {
  W <- as.owin(st_geometry(lek_polygon))
  pts <- st_coordinates(lek_points)
  X <- ppp(pts[,1], pts[,2], window = W)
  
  nn <- nndist(X)
  
  tibble(nnd_mean = mean(nn), nnd_median = median(nn), nnd_sd = sd(nn), nnd_count = length(nn), nnd_cv = sd(nn)/mean(nn))
}

## 2) Effective number of unordered point pairs contributing to each PCF value
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
    
    if (length(w) == 0) return(0)
    
    p <- w / sum(w)
    1 / sum(p^2)
  })
}

## 3) Compute g_inhom(r) per lek x date
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

## 4) Peak detection function
detect_peaks <- function(s, g, med_nnd) {
  
  # Rescale distance to go back to metres
  r <- s * med_nnd
  
  # Apply lower cutoff in rescaled units
  keep <- s >= lower_nnd_bound_s
  s <- s[keep]
  r <- r[keep]
  g <- g[keep]
  
  # Local maxima on the smoothed curve
  is_peak <- g > dplyr::lag(g) & g > dplyr::lead(g)
  peak_idx <- which(is_peak)
  
  # Candidate peaks. Peak height is taken from the smoothed curve.
  peaks <- tibble(idx = peak_idx, s_peak = s[peak_idx], r_peak = r[peak_idx], g_peak = g[peak_idx])
  
  # Prominence window width in rescaled units
  win_half_width <- 0.5
  
  # Step size on rescaled axis for second-derivative curvature
  ds <- median(diff(s), na.rm = TRUE)
  
  # Estimate prominence and curvature from the smoothed curve
  peaks <- purrr::pmap_dfr(peaks, function(idx, s_peak, r_peak, g_peak) {
    
    left_limit_s  <- s_peak - win_half_width
    right_limit_s <- s_peak + win_half_width
    
    left_idx  <- which(s >= left_limit_s & s < s_peak)
    right_idx <- which(s > s_peak & s <= right_limit_s)
    
    # A peak requires support on both sides to define prominence
    left_min  <- min(g[left_idx], na.rm = TRUE)
    right_min <- min(g[right_idx], na.rm = TRUE)
    baseline  <- max(left_min, right_min)
    peak_prominence <- g_peak - baseline
    
    # Curvature: negative second derivative of smoothed g(s)
    gpp <- (g[idx + 1] - 2 * g[idx] + g[idx - 1]) / (ds^2)
    peak_curvature <- -gpp
    
    tibble(idx = idx, s_peak = s_peak, r_peak = r_peak, g_peak = g_peak,
           peak_prominence = peak_prominence, peak_curvature = peak_curvature)
  }) %>%
    filter(!is.na(peak_prominence), peak_prominence >= min_prominence)
  
  if (nrow(peaks) == 0) return(NULL)
  
  # Iteratively enforce minimum peak separation.
  # At each step, retain the most prominent remaining peak and remove
  # all other candidates within min_peak_sep_s on the rescaled distance axis.
  remaining <- peaks %>% arrange(desc(peak_prominence), desc(g_peak))
  selected <- list()
  
  while (nrow(remaining) > 0) {
    chosen <- remaining[1, ]
    selected[[length(selected) + 1]] <- chosen
    
    remaining <- remaining %>% filter(abs(s_peak - chosen$s_peak) >= min_peak_sep_s)
  }
  
  bind_rows(selected) %>% arrange(s_peak) %>% select(-idx)
}

## 5) Compute KDE intensity surface and extract the global mode
compute_intensity_features <- function(lek_polygon, lek_points_sf) {
  
  # Convert points and polygons to spatstat and create a point pattern object
  W <- as.owin(st_geometry(lek_polygon))
  xy <- st_coordinates(lek_points_sf)
  X <- ppp(xy[,1], xy[,2], window = W)
  
  # Estimate KDE bandwidth
  sigma <- bw.ppl(X)
  lambda_hat <- density.ppp(X, sigma = sigma, edge = TRUE, at = "pixels")
  
  # Extract location of maximum intensity (mode)
  v <- lambda_hat$v
  max_v <- max(v, na.rm = TRUE)
  idx_all <- which(v == max_v, arr.ind = TRUE)
  
  idx <- idx_all[1, , drop = FALSE]
  mode <- tibble(mx = lambda_hat$xcol[idx[2]], my = lambda_hat$yrow[idx[1]])
  
  return(mode)
}

## 6) Consider points within the lek polygon only
clip_points_to_polygon <- function(pts_sf, polygon_sf) {
  inside <- lengths(st_within(pts_sf, polygon_sf)) > 0
  pts_sf[inside, ]
}

## 7) Compute a KDE of the points and make a mask that defines the core region
make_kde_grid <- function(pts_sf, lek_polygon, sigma, core_prob = 0.75, dimyx = 256) {
  
  # Convert the polygon to a spatstat window and extract point coordinates
  W <- as.owin(st_geometry(lek_polygon))
  xy <- st_coordinates(pts_sf)
  
  # Build the point pattern and estimate the KDE on a grid
  X <- ppp(x = xy[, 1], y = xy[, 2], window = W)
  dens <- density.ppp(X, sigma = sigma, edge = TRUE, dimyx = dimyx)
  
  # Build the point pattern and estimate the KDE on a grid
  kde_df <- as.data.frame(dens)
  names(kde_df) <- c("x", "y", "z")
  
  # Compute KDE mass and normalize so the KDE values so the mass sums up to 1
  cell_area <- dens$xstep * dens$ystep
  kde_df <- kde_df %>% mutate(p = z / max(z, na.rm = TRUE))
  
  # Mark grid cells as inside or outside the KDE core mask
  kde_df <- kde_df %>% mutate(core_prob = core_prob, sigma_used = sigma, contour_level = 1-core_prob, 
                              in_core = is.finite(p) & (p >= 1-core_prob))
  
  list(kde_df = kde_df, contour_level = 1-core_prob, cell_area = cell_area, xstep = dens$xstep, ystep = dens$ystep)
}

## 8) Check which points fall inside the core region (KDE mask)
subset_points_to_kde_core <- function(pts_sf, kde_grid) {
  
  # Extract point coordinates and only keep points that are within the KDE core
  pts_xy <- st_coordinates(pts_sf)
  core_cells <- kde_grid$kde_df %>% filter(in_core) %>% select(x, y)
  if (nrow(core_cells) == 0) return(pts_sf[0, ])
  
  # Match each point to nearest grid-cell centre in x and y separately
  x_vals <- sort(unique(kde_grid$kde_df$x))
  y_vals <- sort(unique(kde_grid$kde_df$y))
  
  # Get grid cell centres along each axis
  nearest_x <- vapply(pts_xy[, 1], function(xx) x_vals[which.min(abs(x_vals - xx))], numeric(1))
  nearest_y <- vapply(pts_xy[, 2], function(yy) y_vals[which.min(abs(y_vals - yy))], numeric(1))
  
  pts_key <- tibble(x = nearest_x, y = nearest_y)
  core_key <- core_cells %>% mutate(in_core = TRUE)
  
  # Mark and keep points whose matched grid cell lies inside the KDE core
  inside <- pts_key %>% left_join(core_key, by = c("x", "y")) %>%
    mutate(in_core = ifelse(is.na(in_core), FALSE, in_core)) %>% pull(in_core)
  
  pts_sf[inside, ]
}

## 9) Plot current-year point pattern against the previous-year KDE core
plot_crossyear_core_overlap <- function(kde_grid, lek_polygon, pts_prev, pts_curr, pts_curr_in_core,
                                        lek, date_prev, date_curr, plot_dir, plot_limits) {
  
  # Keep only KDE grid cells that define the previous-year core
  core_df <- kde_grid$kde_df %>% filter(in_core)
  
  # Make a safe file name for this transition
  safe_lek <- gsub("[^A-Za-z0-9]+", "_", lek)
  out_file <- file.path(plot_dir, paste0("core_overlap_", safe_lek, "_",
                                         format(as.Date(date_prev), "%Y%m"), "_to_",
                                         format(as.Date(date_curr), "%Y%m"), ".png"))
  
  # Plot the previous-year core, previous points, current points and points used in the computation
  p <- ggplot() +
    geom_tile(data = core_df, aes(x = x, y = y), width = kde_grid$xstep, height = kde_grid$ystep, fill = "#C49102", alpha = 0.4) +
    geom_sf(data = pts_curr, colour = "black", size = 1.2, alpha = 0.6) +
    geom_sf(data = pts_curr_in_core, colour = "#3FA7A3", size = 1.2) +
    coord_sf(xlim = plot_limits$xlim, ylim = plot_limits$ylim, expand = FALSE) + theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12))
  
  ggsave(out_file, p, width = 5, height = 5, dpi = 300)
}

## 10) Compute distance of points from the KDE mode
distance_to_kde_mode <- function(pts_sf, kde_grid, crs_use) {
  mode_cell <- kde_grid$kde_df %>% filter(is.finite(p)) %>% slice_max(order_by = p, n = 1, with_ties = FALSE)
  mode_pt <- st_as_sf(mode_cell, coords = c("x", "y"), crs = crs_use)
  as.numeric(st_distance(pts_sf, mode_pt))
}

## 11) Estimate KDE bandwidth for a given point pattern
get_kde_sigma <- function(pts_sf, lek_polygon) {
  
  # Convert polygon and points to spatstat objects and create a point pattern object
  W <- as.owin(st_geometry(lek_polygon))
  xy <- st_coordinates(pts_sf)
  
  X <- ppp(x = xy[, 1], y = xy[, 2], window = W)
  
  # Estimate the smoothing bandwidth for the KDE
  sigma <- bw.diggle(X, hmax = diameter(W))
  sigma <- as.numeric(sigma)[1]
  
  return(sigma)
}

## 12) Compute nearest-neighbour distance between points
nnd <- function(from_pts, to_pts) {
  dmat <- st_distance(from_pts, to_pts)
  dmin <- apply(dmat, 1, min)
  
  as.numeric(dmin)
}

## 13) Translate points by a fixed distance in a random direction
translate_points_random <- function(pts_sf, shift_dist) {
  
  # Draw a random direction and shift all points by the same offset
  theta <- runif(1, 0, 2 * pi)
  dx <- shift_dist * cos(theta)
  dy <- shift_dist * sin(theta)
  
  pts_xy <- st_coordinates(pts_sf)
  pts_shift <- tibble(x = pts_xy[, 1] + dx, y = pts_xy[, 2] + dy)
  
  st_as_sf(pts_shift, coords = c("x", "y"), crs = st_crs(pts_sf))
}

## 14) Rotate points by a random angle around their centroid
rotate_points_random <- function(pts_sf, angle_max = pi / 6) {
  
  # Draw a random rotation angle and rotate all points around the centroid
  theta <- runif(1, -angle_max, angle_max)
  pts_xy <- st_coordinates(pts_sf)
  ctr <- colMeans(pts_xy)
  pts_ctr <- sweep(pts_xy, 2, ctr, "-")
  
  rot_mat <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, byrow = TRUE)
  pts_rot <- pts_ctr %*% t(rot_mat)
  pts_rot <- sweep(pts_rot, 2, ctr, "+")
  pts_rot <- tibble(x = pts_rot[, 1], y = pts_rot[, 2])
  
  st_as_sf(pts_rot, coords = c("x", "y"), crs = st_crs(pts_sf))
}

## 15) Translate and rotate points
transform_points_random <- function(pts_sf, shift_dist, angle_max = pi / 6) {
  pts_trans <- translate_points_random(pts_sf, shift_dist = shift_dist)
  rotate_points_random(pts_trans, angle_max = angle_max)
}

## 16) Simulate transformed previous-year points and compute nearest-neighbour distances
simulate_transform_crossyear_nnd <- function(prev_pts, curr_pts, n_sims = 999) {
  
  # Compute the shift distance as half of the median nearest-neighbour distance within the previous image
  prev_dmat <- st_distance(prev_pts, prev_pts)
  diag(prev_dmat) <- Inf
  prev_nnd <- apply(prev_dmat, 1, min)
  shift_dist <- median(as.numeric(prev_nnd)) / 2
  
  # Repeat the translate-rotate randomisation and store point-level NNDs from each simulation
  sim_pointwise <- map_dfr(seq_len(n_sims), function(s) {
    prev_transform <- transform_points_random(prev_pts, shift_dist = shift_dist, angle_max = pi / 6)
    transform_nnd <- nnd(curr_pts, prev_transform)
    
    tibble(sim = s, point_id = seq_along(transform_nnd),
           nnd_to_prev = transform_nnd)
  })
  
  # Summarise point-level NNDs from each simulation
  sim_summary <- sim_pointwise %>%
    group_by(sim) %>%
    summarise(mean_nnd_transform_rand = mean(nnd_to_prev),
              median_nnd_transform_rand = median(nnd_to_prev),
              .groups = "drop")
  
  list(summary = sim_summary, pointwise = sim_pointwise)
}