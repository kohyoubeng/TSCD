###############################################################################
# TaR Plot — EGIG.TSCD Model
#
# This script visualizes TaR duration forecasts alongside
# diurnaly adjusted durations for AAPL using the EGIG.TSCD model.
#
# The procedure:
# • Reads TaR forecasts and observed durations from CSV
# • Constructs shaded risk regions for multiple TaR levels:
#     – 90% TaR
#     – 95% TaR
#     – 97.5% TaR
# • Draws boundary lines for each TaR level
# • Plots observed durations as a dotted line
# • Produces a publication-ready plot with customized axes, legend,
#   and formatting
#
# The SAME script applies to both AAPL and TSLA datasets.
#
# --- To use for TSLA ---
# • Replace input file:
#     "AAPL Index, ObsDur, fphi & TaR Forecast_EGIG.TSCD.csv" 
#       → "TSLA Index, ObsDur, fphi & TaR Forecast_EGIG.TSCD.csv"
# • Replace output filename:
#     "AAPL TaR EGIG.TSCD Plot.png" → "TSLA TaR EGIG.TSCD Plot.png"
###############################################################################

rm(list = ls())

TaR = read.csv("AAPL Index, ObsDur, fphi & TaR Forecast_EGIG.TSCD.csv",
               header = TRUE)

library(ggplot2)
lins = c("dotted")
cols = c("gray55", "gray75", "gray90")

# EGIG.TSCD
TSCD = data.frame(
  Index = TaR$Index,
  Obs   = TaR$ObsDur,
  LV90  = TaR$fTaR_0.1,
  LV95  = TaR$fTaR_0.05,
  LV97.5= TaR$fTaR_0.025
)

p = ggplot(data = TSCD, aes(x = Index)) +
  
  # TaR areas
  geom_area(aes(y = LV97.5, color = "97.5TaR", fill = "97.5TaR")) +
  geom_area(aes(y = LV95,   color = "95TaR",   fill = "95TaR")) +
  geom_area(aes(y = LV90,   color = "90TaR",   fill = "90TaR")) +
  
  scale_fill_manual(
    values = c("90TaR" = cols[3], "95TaR" = cols[2], "97.5TaR" = cols[1]),
    labels = c("   90% TaR", "   95% TaR", "97.5% TaR")
  ) +
  
  # TaR lines
  geom_line(aes(y = LV97.5), color = "black", linetype = "solid", size = 0.1) +
  geom_line(aes(y = LV95),   color = "black", linetype = "solid", size = 0.1) +
  geom_line(aes(y = LV90),   color = "black", linetype = "solid", size = 0.1) +
  
  # observed duration line
  geom_line(aes(y = Obs, color = "Obs", linetype = "Obs"), size = 0.2) +
  
  scale_color_manual(
    values = c("97.5TaR" = NA, "95TaR" = NA, "90TaR" = NA, "Obs" = "black")
  ) +
  
  scale_linetype_manual(
    values = c("Obs" = lins[1]),
    labels = c(expression("Observed, " * italic("x")["i"]))
  ) +
  
  theme_classic() +
  
  # X-axis: 500, 1000, 1500, ...
  scale_x_continuous(
    expand = c(0, 5),
    name = "Time Index",
    breaks = seq(0, max(TSCD$Index), by = 500)
  ) +
  
  scale_y_continuous(
    expand = expansion(c(0, 0.02)),
    name = "Diurnally Adjusted Duration",
    breaks = seq(from = 0, to = max(TaR$fTaR_0.025), by = 5)
  ) +
  
  theme(
    legend.title = element_blank(),
    legend.box.background = element_rect(color = "black", size = 0.1, fill = NULL),
    legend.position = c(0.9, 0.88),
    legend.spacing.y = unit(-1.4, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.key.width  = unit(8, "mm"),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    panel.background = element_rect(color = "black", size = 0.3)
  ) +
  
  guides(
    fill = guide_legend(order = 1, reverse = TRUE,
                        override.aes = list(colour = "black")),
    linetype = guide_legend(order = 2, reverse = TRUE),
    color = FALSE
  )

# Save plot
ggsave("AAPL TaR EGIG.TSCD Plot.png", plot = p,
       width = 7, height = 4.5, dpi = 600)
