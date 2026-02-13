###############################################################################
# Conditional vs Empirical PDF/CDF/Survival/Hazard (Burr SCD Model)
#
# This script compares the conditional Burr–SCD model functions with the
# empirical Gamma-kernel estimates of the PDF, CDF, Survival, and Hazard
# functions based on standardized durations.
#
# The same script applies to both AAPL and TSLA datasets.
#
# --- To use for TSLA ---
# • Replace "AAPL In-sample obs, SCD & TSCD.csv" with
#   "TSLA In-sample obs, SCD & TSCD.csv"
# • Replace "AAPL Empirical result.csv" with
#   "TSLA Empirical result.csv"
# • Replace all output filenames beginning with "AAPL" to "TSLA"
# • Update the Burr shape parameter value **a** and **q** to the TSLA estimate
#
# Notes:
# • MUST ensure the input CSV contains columns:
#     x.use  (diurnally adjusted durations)
#     Burr    (fitted conditional mean from Burr–SCD model)
# • The empirical result CSV must contain:
#     x_i_ordered, Emp_pdf, Emp_cdf, Emp_sur, Emp_haz
# • No other changes to the model structure or plotting procedure are required.
###############################################################################

rm(list = ls())

# Required Libraries
library(gsl)      # For gamma functions and Bessel functions (ensure installed)
library(ggplot2)  # For plotting
library(pracma)   # For numerical integration (cumtrapz function)

###############################################################################
# Step 1: Load data and parameters

par.read = c(1.35666081585, 2.18603181599)
a = par.read[1]
q = par.read[2]
b = 1 / (q * beta(1+1/a, q-1/a))

# Load durations and scaling factor phi
x_data = read.csv("Adj dur AAPL In-sample (t)_29 Sep-10 Oct 2025.csv", 
                  header = TRUE)
x = x_data$x.use

info_data = read.csv("AAPL In-sample obs, SCD & TSCD.csv", header = TRUE)
phi = info_data$Burr

# Filter positive durations and corresponding phi
valid_idx = which(x > 0 & phi > 0)
x = x[valid_idx]
phi = phi[valid_idx]

# Standardized durations
eps = x / phi

eps_grid = sort(eps)

###############################################################################
###############################################################################
# Step 2: Load empirical results from CSV and interpolate to eps_grid

G.kernel = read.csv("AAPL Empirical result.csv", header = TRUE)

Emp_pdf_interp = approx(G.kernel$x_i_ordered, G.kernel$Emp_pdf,
                        xout = eps_grid, rule = 2)$y
Emp_cdf_interp = approx(G.kernel$x_i_ordered, G.kernel$Emp_cdf,
                        xout = eps_grid, rule = 2)$y
Emp_survival_interp = approx(G.kernel$x_i_ordered, G.kernel$Emp_sur,
                             xout = eps_grid, rule = 2)$y
Emp_hazard_interp = approx(G.kernel$x_i_ordered, G.kernel$Emp_haz,
                           xout = eps_grid, rule = 2)$y

###############################################################################
# Step 3: Define custom conditional PDF and CDF functions

compute_pdf = function(eps, phi, a, q, b) {
  a * q * (eps)^(a-1) / 
    (b^a * ((1 + b^(-a) * (eps)^a))^(1+q)) / phi
}

compute_cdf = function(eps, a, q, b) {
  1 - (1 + (eps * q * beta(1+1/a, q-1/a))^a)^(-q)
}


###############################################################################
# Step 4: Compute conditional model values on grid

# Use median phi for smooth conditional curve
median_phi = median(phi, na.rm = TRUE)

# Use smoother grid to avoid overly smooth plots
eps_gridCond = seq(min(eps, na.rm = TRUE),
                   max(eps, na.rm = TRUE),
                   length.out = length(x))

# Compute conditional PDF, CDF, survival, hazard on the grid
cond_pdf = compute_pdf(eps_gridCond, median_phi, a, q, b)
cond_cdf = compute_cdf(eps_gridCond, a, q, b)
cond_survival = 1 - cond_cdf
cond_hazard = cond_pdf / cond_survival

# Map conditional epsilon grid back to raw scale
x_gridCond = eps_gridCond * median_phi

###############################################################################
# Step 5: Combine empirical and conditional results for plotting

# Interpolate conditional functions on empirical grid for alignment
cond_pdf_interp = approx(x_gridCond, cond_pdf, xout = eps_grid, rule = 2)$y
cond_cdf_interp = approx(x_gridCond, cond_cdf, xout = eps_grid, rule = 2)$y
cond_survival_interp = approx(x_gridCond, cond_survival, xout = eps_grid, 
                              rule = 2)$y
cond_hazard_interp = approx(x_gridCond, cond_hazard, xout = eps_grid, rule = 2)$y

result_df = data.frame(
  eps_i = eps_grid,
  conditional_pdf = cond_pdf_interp,
  conditional_cdf = cond_cdf_interp,
  conditional_survival = cond_survival_interp,
  conditional_hazard = cond_hazard_interp,
  empirical_pdf = Emp_pdf_interp,
  empirical_cdf = Emp_cdf_interp,
  empirical_survival = Emp_survival_interp,
  empirical_hazard = Emp_hazard_interp
)

write.csv(result_df, "AAPL Cond vs Emp_Burr.SCD.csv", row.names = FALSE)

###############################################################################
# Step: Truncate or extend data to x = 20

extend_to_cutoff = function(df, cutoff = 20) {
  y_cols = setdiff(names(df), "eps_i")
  max_eps = max(df$eps_i, na.rm = TRUE)
  
  if (max_eps > cutoff) {
    df_trim = df[df$eps_i <= cutoff, ]
    # Interpolate each column to get value at cutoff
    last_values = sapply(y_cols, function(col) {
      approx(df$eps_i, df[[col]], xout = cutoff, rule = 2)$y
    })
    last_row = data.frame(eps_i = cutoff, as.list(last_values))
    df_trim = rbind(df_trim, last_row)
  } else {
    df_trim = df
    # Extend by repeating last row at cutoff
    last_row = df_trim[nrow(df_trim), ]
    last_row$eps_i = cutoff
    df_trim = rbind(df_trim, last_row)
  }
  
  return(df_trim)
}

# Apply it to result_df before plotting
filtered_data = extend_to_cutoff(result_df, cutoff = 20)

###############################################################################
# Step 6: Plotting

common_theme = theme_minimal(base_size = 15) +
  theme(
    legend.title = element_blank(),
    legend.box.background = element_rect(color = "black", size = 0.1, 
                                         linetype = "solid"),
    legend.spacing.y = unit(-1.4, "mm"),
    legend.key.height = unit(3, 'mm'),
    legend.key.width = unit(10, 'mm'),
    legend.text = element_text(size = 16),
    legend.background = element_blank(),
    panel.background = element_rect(color = NULL, size = 0.3),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text = element_text(size = 20),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

x_scale = scale_x_continuous(breaks = seq(0, 20, by = 1), limits = c(0, 20))
y_scale = scale_y_continuous(breaks = seq(0, 1, by = 0.2))

line_type_settings = scale_linetype_manual(
  values = c("Empirical" = "solid", "SCD_Burr model" = "dotted"),
  labels = c("Empirical", expression("SCD"[Burr] * " model"))
)

# PDF plot
pdf_plot = ggplot(filtered_data, aes(x = eps_i)) +
  geom_line(aes(y = empirical_pdf, linetype = "Empirical"), color="black") +
  geom_line(aes(y = conditional_pdf, linetype = "SCD_Burr model"), color="blue") +
  labs(title = "PDF", x = expression(epsilon[i]), y = "Pdf") +
  line_type_settings + common_theme + x_scale + 
  y_scale +
  coord_cartesian(ylim = c(0, max(c(filtered_data$empirical_pdf, filtered_data$conditional_pdf), na.rm = TRUE) * 1.1)) +
  theme(legend.position = c(0.88, 0.73))

# CDF plot
cdf_plot = ggplot(filtered_data, aes(x = eps_i)) +
  geom_line(aes(y = empirical_cdf, linetype = "Empirical"), color="black") +
  geom_line(aes(y = conditional_cdf, linetype = "SCD_Burr model"), color="blue") +
  labs(title = "CDF", x = expression(epsilon[i]), y = "Cdf") +
  line_type_settings + common_theme + x_scale + 
  y_scale +
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.position = c(0.88, 0.60))

# Survival plot
survival_plot = ggplot(filtered_data, aes(x = eps_i)) +
  geom_line(aes(y = empirical_survival, linetype = "Empirical"), color="black") +
  geom_line(aes(y = conditional_survival, linetype = "SCD_Burr model"), color="blue") +
  labs(title = "Survival", x = expression(epsilon[i]), y = "Survival") +
  line_type_settings + common_theme + x_scale + y_scale +
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.position = c(0.88, 0.73))

# Hazard plot
hazard_plot = ggplot(filtered_data, aes(x = eps_i)) +
  geom_line(aes(y = empirical_hazard, linetype = "Empirical"), color="black") +
  geom_line(aes(y = conditional_hazard, linetype = "SCD_Burr model"), color="blue") +
  labs(title = "Hazard", x = expression(epsilon[i]), y = "Hazard") +
  line_type_settings + common_theme + x_scale + 
  y_scale +
  coord_cartesian(ylim = c(0, max(c(filtered_data$empirical_hazard, 
                                    filtered_data$conditional_hazard), 
                                  na.rm = TRUE) * 1.1)) +
  theme(legend.position = c(0.88, 0.73))

# Show plots
pdf_plot
cdf_plot
survival_plot
hazard_plot

# Save combined plots
library(gridExtra)
combined_plot = grid.arrange(pdf_plot, cdf_plot, survival_plot, hazard_plot, 
                             ncol = 1)
ggsave("AAPL combined Cond plot_Burr.SCD.png", plot = combined_plot, width = 16, 
       height = 12, dpi = 300)

# Save hazard plot individually with no title and smaller font
hazard_plot_no_title = hazard_plot + 
  labs(title = NULL) +
  theme(
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )

ggsave("AAPL Hazard_Burr.SCD.png", plot = hazard_plot_no_title, width = 16, 
       height = 3, dpi = 300)
