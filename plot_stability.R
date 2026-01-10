rm(list = ls())

## Load libraries
library(tidyverse)

## Save figure parameters
save_pub_fig <- function(plot, filename_base, width = 7, height = 5) {
  ggsave(paste0(filename_base, ".png"),
         plot = plot,
         width = width,
         height = height,
         units = "in",
         dpi = 600)
}

## Load stability metrics
stab <- read.csv("derived_metrics/stability_ALL.csv")

## Factor ordering + labels
stab$lek_id <- factor(stab$lek_id,
                      levels = c("Velavadar_LEK1",
                                 "Velavadar_LEK2",
                                 "TalChhapar_TC"),
                      labels = c("Velavadar Lek 1",
                                 "Velavadar Lek 2",
                                 "Tal Chhapar"))

## Colour palettes
fill_cols <- c(
  "Velavadar Lek 1" = "#4DAF4A",
  "Velavadar Lek 2" = "#377EB8",
  "Tal Chhapar"     = "#D6604D"
)

point_cols <- c(
  "Velavadar Lek 1" = "#1B7837",
  "Velavadar Lek 2" = "#2166AC",
  "Tal Chhapar"     = "#8B1A1A"
)

## Centroid displacement
p_centroid <- ggplot(stab, aes(x = lek_id, y = centroid_shift)) +
  geom_boxplot(aes(fill = lek_id),
               width = 0.55,
               outlier.shape = NA,
               linewidth = 0.8,
               colour = "black") +
  geom_jitter(aes(colour = lek_id),
              width = 0.10,
              size = 2.4,
              alpha = 0.9) +
  scale_fill_manual(values = fill_cols) +
  scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Lek", y = "Intensity centroid displacement (m)")

save_pub_fig(p_centroid, "fig_stability_centroid_shift")

## Mode displacement
p_mode <- ggplot(stab, aes(x = lek_id, y = mode_shift)) +
  geom_boxplot(aes(fill = lek_id),
               width = 0.55,
               outlier.shape = NA,
               linewidth = 0.8,
               colour = "black") +
  geom_jitter(aes(colour = lek_id),
              width = 0.10,
              size = 2.4,
              alpha = 0.9) +
  scale_fill_manual(values = fill_cols) +
  scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Lek", y = "Intensity mode displacement (m)")

save_pub_fig(p_mode, "fig_stability_mode_shift")

## Cross-year NN distance
p_cross_nn <- ggplot(stab, aes(x = lek_id, y = nn_cross_median)) +
  geom_boxplot(aes(fill = lek_id),
               width = 0.55,
               outlier.shape = NA,
               linewidth = 0.8,
               colour = "black") +
  geom_jitter(aes(colour = lek_id),
              width = 0.10,
              size = 2.4,
              alpha = 0.9) +
  scale_fill_manual(values = fill_cols) +
  scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Lek", y = "Cross-year nearest-neighbour distance (m)")

save_pub_fig(p_cross_nn, "fig_stability_crossyear_nn")

## ---- Time series: cross-year NN distance ----

p_nn_ts <- ggplot(stab, aes(x = date_now, y = nn_cross_median, colour = lek_id, group = lek_id)) +
  geom_line(linewidth = 1) + geom_point(size = 2.6) + scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Year", y = "Cross-year nearest-neighbour distance (m)", colour = "Lek")

save_pub_fig(p_nn_ts, "fig_stability_crossyear_nn_timeseries", width = 8, height = 4.8)

## ---- Time series: centroid displacement ----

p_centroid_ts <- ggplot(stab, aes(x = date_now, y = centroid_shift, colour = lek_id, group = lek_id)) +
  geom_line(linewidth = 1) + geom_point(size = 2.6) + scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Year", y = "Intensity centroid displacement (m)", colour = "Lek")

save_pub_fig(p_centroid_ts, "fig_stability_centroid_timeseries", width = 8, height = 4.8)

## ---- Time series: mode displacement ----

p_modes_ts <- ggplot(stab, aes(x = date_now, y = mode_shift, colour = lek_id, group = lek_id)) +
  geom_line(linewidth = 1) + geom_point(size = 2.6) + scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Year", y = "Intensity centroid displacement (m)", colour = "Lek")

save_pub_fig(p_modes_ts, "fig_stability_mode_timeseries", width = 8, height = 4.8)
