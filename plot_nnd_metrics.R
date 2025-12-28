rm(list = ls())

## Load libraries
library(tidyverse)

## Save figure parameters
save_pub_fig <- function(plot, filename_base, width = 7, height = 5) {
  ggsave(paste0(filename_base, ".png"), plot = plot, 
         width = width, height = height, units = "in", dpi = 600)
}

## Load NND metrics (NEW workflow)
nnd_all <- read.csv("derived_metrics/nnd_ALL.csv")

## Factor ordering + labels
nnd_all$lek_id <- factor(nnd_all$lek_id,
  levels = c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC"),
  labels = c("Velavadar Lek 1", "Velavadar Lek 2", "Tal Chhapar"))

## Colour palettes
fill_cols <- c(
  "Velavadar Lek 1" = "#4DAF4A",  # manuscript green
  "Velavadar Lek 2" = "#377EB8",  # manuscript blue
  "Tal Chhapar"     = "#D6604D"   # muted red fill
)

point_cols <- c(
  "Velavadar Lek 1" = "#1B7837",  # dark green
  "Velavadar Lek 2" = "#2166AC",  # dark blue
  "Tal Chhapar"     = "#8B1A1A"   # dark red
)

# Mean NND
p_mean <- ggplot(nnd_all, aes(x = lek_id, y = nnd_mean)) +
  geom_boxplot(aes(fill = lek_id), width = 0.55, outlier.shape = NA, linewidth = 0.8, colour = "black") +
  geom_jitter(aes(colour = lek_id), width = 0.10, size = 2.4, alpha = 0.9) +
  scale_fill_manual(values = fill_cols) + scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Lek", y = "Mean nearest-neighbour distance (m)")

save_pub_fig(p_mean, "fig_nnd_mean_by_lek")

# SD of NND (key regularity metric)
p_sd <- ggplot(nnd_all, aes(x = lek_id, y = nnd_sd)) +
  geom_boxplot(aes(fill = lek_id), width = 0.55, outlier.shape = NA, linewidth = 0.8, colour = "black") +
  geom_jitter(aes(colour = lek_id), width = 0.10, size = 2.4, alpha = 0.9) +
  scale_fill_manual(values = fill_cols) + scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", 
        axis.title.x = element_text(margin = margin(t = 10)), 
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Lek", y = "Standard deviation of nearest-neighbour distance (m)")

save_pub_fig(p_sd, "fig_nnd_sd_by_lek")

# NND count
p_count <- ggplot(nnd_all, aes(x = lek_id, y = nnd_count)) +
  geom_boxplot(aes(fill = lek_id), width = 0.55, outlier.shape = NA, linewidth = 0.8, colour = "black") +
  geom_jitter(aes(colour = lek_id), width = 0.10, size = 2.4, alpha = 0.9) +
  scale_fill_manual(values = fill_cols) + scale_colour_manual(values = point_cols) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Lek", y = "Number of territories")

save_pub_fig(p_count, "fig_nnd_count_by_lek")
