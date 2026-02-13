###############################################################################
# Conditional vs Empirical PDF/CDF/Survival/Hazard (EGIG TSCD Model)
#
# This script compares the conditional EGIG–TSCD model functions with the
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
# • Update the EGIG shape parameter value **lam1**, **lam2**, **del1**,  
#   **del2**, **w1** and **w2** to the TSLA estimate
#
# Notes:
# • MUST ensure the input CSV contains columns:
#     x.use  (diurnally adjusted durations)
#     TSCD    (fitted conditional mean from EGIG–TSCD model)
# • The empirical result CSV must contain:
#     x_i_ordered, Emp_pdf, Emp_cdf, Emp_sur, Emp_haz
# • No other changes to the model structure or plotting procedure are required.
###############################################################################

rm(list = ls())

# Required Libraries
library(gsl)
library(ggplot2)
library(pracma)
library(gridExtra)

###############################################################################
# Step 1: Load data and two-regime parameters (lam1, lam2, del1, del2, w1, w2)
par.read = c(0.054160305264, 1.075434255652, 
             1.820562901656, 0.457217598548, 
             0.324900521074, 0.093415855822)

lam1 = par.read[1]; lam2 = par.read[2]
del1 = par.read[3]; del2 = par.read[4]
w1   = par.read[5]; w2   = par.read[6]

# Pre-compute Bessel values for both regimes
k1_regime1 = besselK(nu = (lam1 + 1) / del1, x = w1)
k2_regime1 = besselK(nu = lam1 / del1,       x = w1)
k1_regime2 = besselK(nu = (lam2 + 1) / del2, x = w2)
k2_regime2 = besselK(nu = lam2 / del2,       x = w2)

info_data = read.csv("AAPL In-sample obs, SCD & TSCD.csv", header = TRUE)
x = info_data$x.use
phi = info_data$TSCD

# standardized durations per-observation
eps = x / phi 

eps_grid = sort(eps)

# threshold (as you requested)
r = median(x, na.rm = TRUE)   # threshold in raw x-space: Regime1 if x <= r, Regime2 if x > r

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
# Step 3: Define conditional PDF and CDF (per-observation) — robust integrate with tryCatch

compute_pdf = function(eps, phi_i, lam, de, w, k1, k2) {
  # eps here is standardized (x_i/phi_i)
  ((de * eps^(lam - 1) * k1^lam) / (2 * k2^(lam + 1)) *
     exp(-w / 2 * ((eps * k1 / k2)^de + (eps * k1 / k2)^(-de)))) / phi_i
}

generalized_gamma_inc_safe = function(eta, z, c) {
  # z can be vector; return same length
  sapply(z, function(zi) {
    if (!is.finite(zi) || zi <= 0) return(0)
    val = tryCatch({
      integrate(function(t) t^(eta - 1) * exp(-t - c / t),
                lower = 0, upper = zi, rel.tol = 1e-8, abs.tol = 0)$value
    }, error = function(e) NA_real_)
    if (is.na(val)) {
      # fallback small approximation
      return(NA_real_)
    } else return(val)
  })
}

compute_cdf = function(eps, lam, de, w, k1, k2) {
  # returns cdf on standardized scale
  eta = lam / de
  z = (w / 2) * (eps * k1 / k2)^de
  c = w^2 / 4
  gamma_gen_value = generalized_gamma_inc_safe(eta, z, c)
  denom = 2^(1 - eta) * w^eta * k2
  cdf_val = gamma_gen_value / denom
  # clamp to [0,1] in case of small numeric issues
  cdf_val[!is.finite(cdf_val)] = NA_real_
  cdf_val = pmin(pmax(cdf_val, 0), 1)
  return(cdf_val)
}

###############################################################################
# Step 4 (Revised): Smooth conditional model over a regular eps grid

# Smooth grid for standardized epsilon
eps_gridCond = seq(min(x/phi, na.rm = TRUE),
                   max(x/phi, na.rm = TRUE),
                   length.out = length(x))

# Map to raw x scale using median phi (for consistent scaling)
median_phi = median(phi, na.rm = TRUE)
x_gridCond = eps_gridCond * median_phi

# Regime assignment based on x_gridCond <= r
is_r1 = x_gridCond <= r
is_r2 = x_gridCond > r

# Compute conditional model values (smooth, regime-based)
cond_pdf_r1 = compute_pdf(eps_gridCond, median_phi, lam1, del1, w1, k1_regime1, 
                          k2_regime1)
cond_cdf_r1 = compute_cdf(eps_gridCond, lam1, del1, w1, k1_regime1, k2_regime1)
cond_sur_r1 = 1 - cond_cdf_r1
cond_haz_r1 = cond_pdf_r1 / cond_sur_r1
cond_haz_r1[!is.finite(cond_haz_r1)] = NA

cond_pdf_r2 = compute_pdf(eps_gridCond, median_phi, lam2, del2, w2, k1_regime2, 
                          k2_regime2)
cond_cdf_r2 = compute_cdf(eps_gridCond, lam2, del2, w2, k1_regime2, k2_regime2)
cond_sur_r2 = 1 - cond_cdf_r2
cond_haz_r2 = cond_pdf_r2 / cond_sur_r2
cond_haz_r2[!is.finite(cond_haz_r2)] = NA

# Merge both regimes for plotting
cond_pdf_i = ifelse(is_r1, cond_pdf_r1, cond_pdf_r2)
cond_cdf_i = ifelse(is_r1, cond_cdf_r1, cond_cdf_r2)
cond_sur_i = ifelse(is_r1, cond_sur_r1, cond_sur_r2)
cond_haz_i = ifelse(is_r1, cond_haz_r1, cond_haz_r2)

###############################################################################
# Step 5: Align conditional model with empirical eps_grid via interpolation
cond_pdf_i = approx(x = eps_gridCond, y = cond_pdf_i, xout = eps_grid, rule = 2)$y
cond_cdf_i = approx(x = eps_gridCond, y = cond_cdf_i, xout = eps_grid, rule = 2)$y
cond_sur_i = approx(x = eps_gridCond, y = cond_sur_i, xout = eps_grid, rule = 2)$y
cond_haz_i = approx(x = eps_gridCond, y = cond_haz_i, xout = eps_grid, rule = 2)$y

# Split for color plotting
eps_r = median_phi * eps_gridCond[which.min(abs(x_gridCond - r))]
is_r1_grid = eps_grid <= eps_r
is_r2_grid = eps_grid > eps_r

cond_pdf_r1_plot = ifelse(is_r1_grid, cond_pdf_i, NA)
cond_pdf_r2_plot = ifelse(is_r2_grid, cond_pdf_i, NA)
cond_cdf_r1_plot = ifelse(is_r1_grid, cond_cdf_i, NA)
cond_cdf_r2_plot = ifelse(is_r2_grid, cond_cdf_i, NA)
cond_sur_r1_plot = ifelse(is_r1_grid, cond_sur_i, NA)
cond_sur_r2_plot = ifelse(is_r2_grid, cond_sur_i, NA)
cond_haz_r1_plot = ifelse(is_r1_grid, cond_haz_i, NA)
cond_haz_r2_plot = ifelse(is_r2_grid, cond_haz_i, NA)


###############################################################################
# Step 6: Save result CSV (same name you expected)
result_df = data.frame(
  eps_i = eps_grid,
  conditional_pdf_r1 = cond_pdf_r1_plot,
  conditional_pdf_r2 = cond_pdf_r2_plot,
  conditional_cdf_r1 = cond_cdf_r1_plot,
  conditional_cdf_r2 = cond_cdf_r2_plot,
  conditional_survival_r1 = cond_sur_r1_plot,
  conditional_survival_r2 = cond_sur_r2_plot,
  conditional_hazard_r1 = cond_haz_r1_plot,
  conditional_hazard_r2 = cond_haz_r2_plot,
  empirical_pdf = Emp_pdf_interp,
  empirical_cdf = Emp_cdf_interp,
  empirical_survival = Emp_survival_interp,
  empirical_hazard = Emp_hazard_interp
)

write.csv(result_df, "AAPL Cond vs Emp_EGIG.TSCD.csv", row.names = FALSE)

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
# Step 7: Plotting (kept your formatting)
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
  values = c("Empirical" = "solid",
             "TSCD_EGIG Regime 1" = "dotted",
             "TSCD_EGIG Regime 2" = "dotdash"),
  labels = c("Empirical",
             expression("TSCD"[EGIG] * " Regime 1"),
             expression("TSCD"[EGIG] * " Regime 2"))
)


# PDF
pdf_plot = ggplot(filtered_data, aes(x = eps_i)) +
  geom_line(aes(y = empirical_pdf, linetype = "Empirical"), color = "black", 
            size = 1.0) +
  geom_line(aes(y = conditional_pdf_r1, linetype = "TSCD_EGIG Regime 1"), 
            color = "blue", size = 1.0, na.rm = TRUE) +
  geom_line(aes(y = conditional_pdf_r2, linetype = "TSCD_EGIG Regime 2"), 
            color = "red",  size = 1.0, na.rm = TRUE) +
  labs(title = "PDF", x = expression(epsilon[i]), y = "Pdf") +
  line_type_settings + common_theme + x_scale + 
  y_scale +
  coord_cartesian(ylim = c(0, max(filtered_data$empirical_pdf, na.rm = TRUE, finite = TRUE) * 1.1)) +
  theme(legend.position = c(0.86, 0.73))

# CDF
cdf_plot = ggplot(filtered_data, aes(x = eps_i)) +
  geom_line(aes(y = empirical_cdf, linetype = "Empirical"), color = "black", 
            size = 1.0) +
  geom_line(aes(y = conditional_cdf_r1, linetype = "TSCD_EGIG Regime 1"), 
            color = "blue", size = 1.0, na.rm = TRUE) +
  geom_line(aes(y = conditional_cdf_r2, linetype = "TSCD_EGIG Regime 2"), 
            color = "red",  size = 1.0, na.rm = TRUE) +
  labs(title = "CDF", x = expression(epsilon[i]), y = "Cdf") +
  line_type_settings + common_theme + x_scale + 
  y_scale +
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.position = c(0.86, 0.60))

# Survival
survival_plot = ggplot(filtered_data, aes(x = eps_i)) +
  geom_line(aes(y = empirical_survival, linetype = "Empirical"), 
            color = "black", size = 1.0) +
  geom_line(aes(y = conditional_survival_r1, linetype = "TSCD_EGIG Regime 1"), 
            color = "blue", size = 1.0, na.rm = TRUE) +
  geom_line(aes(y = conditional_survival_r2, linetype = "TSCD_EGIG Regime 2"), 
            color = "red",  size = 1.0, na.rm = TRUE) +
  labs(title = "Survival", x = expression(epsilon[i]), y = "Survival") +
  line_type_settings + common_theme + x_scale + 
  y_scale +
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.position = c(0.86, 0.73))

# Hazard
hazard_plot = ggplot(filtered_data, aes(x = eps_i)) +
  geom_line(aes(y = empirical_hazard, linetype = "Empirical"), color = "black",
            size = 1.0) +
  geom_line(aes(y = conditional_hazard_r1, linetype = "TSCD_EGIG Regime 1"), 
            color = "blue", size = 1.0, na.rm = TRUE) +
  geom_line(aes(y = conditional_hazard_r2, linetype = "TSCD_EGIG Regime 2"), 
            color = "red",  size = 1.0, na.rm = TRUE) +
  labs(title = "Hazard", x = expression(epsilon[i]), y = "Hazard") +
  line_type_settings + common_theme + x_scale + 
  y_scale +
  coord_cartesian(ylim = c(0, max(filtered_data$empirical_hazard, na.rm = TRUE,
                                  finite = TRUE) * 1.1)) +
  theme(legend.position = c(0.86, 0.73))

# Display and save
combined_plot = grid.arrange(pdf_plot, cdf_plot, survival_plot, hazard_plot,
                             ncol = 1)
ggsave("AAPL combined Cond plot_EGIG.TSCD.png", plot = combined_plot, 
       width = 16, height = 12, dpi = 300)

# Save hazard-only (no title)
hazard_plot_no_title = hazard_plot + 
  labs(title = NULL) +
  theme(
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )

ggsave("AAPL Hazard_EGIG.TSCD.png", plot = hazard_plot_no_title, width = 16, 
       height = 3, dpi = 300)

