rm(list = ls())

## Load libraries
library(tidyverse)

## Save figure parameters
save_pub_fig <- function(plot, filename_base, width = 7, height = 5) {
  ggsave(paste0(filename_base, ".png"), plot = plot,
         width = width, height = height, units = "in", dpi = 600)
}

## Load PCF curves
pcf_all <- read.csv("derived_metrics/pcf_curve_ALL.csv")

## Factor ordering + labels (identical to NND plots)
pcf_all$lek_id <- factor(pcf_all$lek_id,
  levels = c("Velavadar_LEK1", "Velavadar_LEK2", "TalChhapar_TC"),
  labels = c("Velavadar Lek 1", "Velavadar Lek 2", "Tal Chhapar"))

## Colour palettes
fill_cols <- c(
  "Velavadar Lek 1" = "#4DAF4A",  # manuscript green
  "Velavadar Lek 2" = "#377EB8",  # manuscript blue
  "Tal Chhapar"     = "#D6604D"   # muted red fill
)

line_cols <- c(
  "Velavadar Lek 1" = "#1B7837",  # dark green
  "Velavadar Lek 2" = "#2166AC",  # dark blue
  "Tal Chhapar"     = "#8B1A1A"   # dark red
)

## Summarise PCF curves across dates within each lek
pcf_summary <- pcf_all %>% group_by(lek_id, r) %>% summarise(
    g_mean = mean(g, na.rm = TRUE), g_sd = sd(g, na.rm = TRUE),
    g_lo = g_mean - g_sd, g_hi = g_mean + g_sd, .groups = "drop")

## Plot mean g(r) with ribbons
p_pcf <- ggplot(pcf_summary, aes(x = r, y = g_mean)) +
  ## Ribbon shows temporal variability within a lek
  geom_ribbon(aes(ymin = g_lo, ymax = g_hi, fill = lek_id), alpha = 0.25, linewidth = 0) +
  geom_line(aes(colour = lek_id), linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, colour = "black") +
  scale_fill_manual(values = fill_cols) + scale_colour_manual(values = line_cols) +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Distance r (m)", y = expression(g[inhom](r)))

## Save figure
save_pub_fig(p_pcf, "fig_pcf_mean_with_ribbon")
