###############################################################################
# Empirical Functions & Plots — Gamma Kernel (Fully Empirical)
#
# This script computes the empirical PDF, CDF, Survival, and Hazard functions
# using Gamma kernel density estimation and generates the corresponding plots.
#
# The same script applies to both AAPL and TSLA datasets.
#
# To use for TSLA:
#   • Replace "AAPL In-sample obs, SCD & TSCD.csv" with
#     "TSLA In-sample obs, SCD & TSCD.csv"
#   • Replace all output filenames beginning with "AAPL" to "TSLA"
#
# Notes:
#   • MUST ensure the input CSV file contains the column "x.use"
#     (diurnally adjusted durations).
#   • The script structure and estimation procedure remain unchanged.
###############################################################################

rm(list = ls())

library(ggplot2)
library(pracma)
library(gridExtra)

# --- Load data (empirical durations only) ---
x_i_data1 = read.csv("AAPL In-sample obs, SCD & TSCD.csv", header = TRUE)
x_i_ordered = sort(x_i_data1$x.use)

###############################################################################
# Bandwidth & Grid
bandwidth = ((0.9 * sd(x_i_ordered) * length(x_i_ordered))^(-0.2))^2
grid_points = seq(0, max(x_i_ordered), length.out = length(x_i_ordered))

# --- Gamma kernel density estimation function ---
gamma_kde = function(data, grid_points, bandwidth) {
  n = length(data)
  Est_pdf = numeric(length(grid_points))
  for (i in seq_along(grid_points)) {
    kernel_values = dgamma(grid_points[i] / data, shape = 1 / bandwidth, scale = bandwidth)
    Est_pdf[i] = mean(kernel_values / data)
  }
  return(Est_pdf)
}

# --- Compute empirical PDF, CDF, Survival, Hazard ---
Emp_pdf = gamma_kde(x_i_ordered, grid_points, bandwidth)
Emp_cdf = pracma::cumtrapz(grid_points, Emp_pdf)
if (tail(Emp_cdf, 1) > 1) Emp_cdf = Emp_cdf / max(Emp_cdf)
Emp_sur = 1 - Emp_cdf
Emp_haz = ifelse(Emp_sur > 0, Emp_pdf / Emp_sur, 0)

# --- Combine and save results ---
results = data.frame(
  x_i_ordered = grid_points,
  Emp_pdf = Emp_pdf,
  Emp_cdf = Emp_cdf,
  Emp_sur = Emp_sur,
  Emp_haz = Emp_haz
)
write.csv(results, "AAPL  Empirical result.csv", row.names = FALSE)

###############################################################################
# Plot Graphs
###############################################################################
common_theme <- theme_minimal(base_size = 15) +
  theme(
    legend.title = element_blank(),
    legend.box.background = element_rect(color = "black", size = 0.1, linetype = "solid"),
    legend.spacing.y = unit(-1.4, "mm"),
    legend.key.height = unit(3, 'mm'),
    legend.key.width = unit(8, 'mm'),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    panel.background = element_rect(color = "black", size = 0.3),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.3, "cm"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

# --- Axes settings ---
x_breaks = seq(0, 56, by = 2)  # 2,4,6,...,56
x_scale = scale_x_continuous(breaks = x_breaks, limits = c(0, max(results$x_i_ordered)))

# Compute ymax to avoid clipping; ensure tick sequence up to 1 is possible
ymax_pdf = ceiling(max(results$Emp_pdf, na.rm = TRUE) * 10) / 10
ymax_haz = ceiling(max(results$Emp_haz, na.rm = TRUE) * 10) / 10

# ensure at least 1 is included in tick sequence
ymax_pdf_plot = max(ymax_pdf, 1)
ymax_haz_plot = max(ymax_haz, 1)

###############################################################################
# --- Plot PDF ---
Emp.pdf.plot = ggplot(results) +
  geom_line(aes(x = x_i_ordered, y = Emp_pdf), color = "black", na.rm = TRUE) +
  labs(
    title = "Empirical PDF (Gamma-KDE)",
    x = "x",
    y = "Empirical PDF"
  ) +
  common_theme + x_scale +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),            # numbering labels 0.2,0.4,...,1.0
    limits = c(0, ymax_pdf_plot)            # include 1 but don't clip curve
  )

# --- Plot CDF ---
Emp.cdf.plot = ggplot(results) +
  geom_line(aes(x = x_i_ordered, y = Emp_cdf), color = "black", na.rm = TRUE) +
  labs(
    title = "Empirical CDF (Gamma-KDE)",
    x = "x",
    y = "Empirical CDF"
  ) +
  common_theme + x_scale +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1)
  )

# --- Plot Survival ---
Emp.sur.plot = ggplot(results) +
  geom_line(aes(x = x_i_ordered, y = Emp_sur), color = "black", na.rm = TRUE) +
  labs(
    title = "Empirical Survival (Gamma-KDE)",
    x = "x",
    y = "Empirical Survival"
  ) +
  common_theme + x_scale +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1)
  )

# --- Plot Hazard ---
Emp.haz.plot = ggplot(results) +
  geom_line(aes(x = x_i_ordered, y = Emp_haz), color = "black", na.rm = TRUE) +
  labs(
    title = "Empirical Hazard (Gamma-KDE)",
    x = "x",
    y = "Empirical Hazard"
  ) +
  common_theme + x_scale +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),            # numbering labels 0.2,0.4,...,1.0
    limits = c(0, ymax_haz_plot)            # include 1 but don't clip curve
  )

###############################################################################
# Display & Save Combined Plot
combined_Emp_plot = gridExtra::grid.arrange(
  Emp.pdf.plot, Emp.cdf.plot, Emp.sur.plot, Emp.haz.plot, ncol = 1
)

ggsave("AAPL  Combined empirical plot.png",
       plot = combined_Emp_plot, width = 16, height = 12, dpi = 300)

###############################################################################
# Save Single Hazard Plot
Emp.haz.plot_no_title = Emp.haz.plot +
  labs(title = NULL) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

ggsave("AAPL Hazard empirical.png",
       plot = Emp.haz.plot_no_title, width = 10, height = 3, dpi = 300)
