###############################################################################
# Plot Observed vs Fitted Durations for Multiple Models
#
# This script generates plots of observed durations and fitted SCD/TSCD models:
#   • Reads in the in-sample data for a chosen stock (AAPL or TSLA)
#   • Plots observed durations (grey line) and fitted model durations (black for individual models)
#   • Saves individual model plots as "STOCK_obs_vs_<ModelName>_fitted.png"
#   • Plots all fitted models together with distinct colors and line styles
#   • Saves combined plot as "STOCK_obs_and_all_models_fitted.png"
#
# --- To apply to AAPL or TSLA ---
# • Change the input CSV filename in the read.csv() line to the desired stock's data
# • All other code works identically for both datasets
#
# Notes:
# • MUST get all models' fitted values for AAPL or TSLA and save in
#    "AAPL In-sample obs, SCD & TSCD.csv" or "TSLA In-sample obs, SCD & TSCD.csv"
#    with column headings "x.use, Wei, GG, Burr, GB2, EGIG, TSCD"
# • The script automatically adjusts the Y-axis to the data range for better visualization
# • Model names, colors, and line styles are predefined but can be modified if needed
###############################################################################

# --- Clear environment ---
rm(list = ls())

# --- Read data ---
data = read.csv("AAPL In-sample obs, SCD & TSCD.csv", header = TRUE)
n = nrow(data)

# --- Define models and styles ---
models = c("Wei", "GG", "Burr", "GB2", "EGIG", "TSCD")
cols_all = c("red", "blue", "darkgreen", "purple", "orange", "black")
ltys_all = c("dotted", "dotted", "dotted", "dotted", "dashed", "solid")

# --- Common axis setup ---
x_axis_seq = seq(0, n, 500)
y_max = max(data[, c("x.use", models)], na.rm = TRUE)  # <-- get true max across all columns
y_axis_seq = pretty(c(0, y_max), n = 11)                # <-- use pretty() for nice, consistent ticks

# --- 1. Plot and save individual models (grey obs + black fitted) ---
for (model_name in models) {
  png(filename = paste0("AAPL_obs_vs_", model_name, "_fitted.png"),
      width = 1200, height = 500, res = 150)
  
  par(mar = c(4, 4, 1, 1))
  
  # Plot observed durations
  plot(data$x.use, type = "l", col = "grey", lty = "solid",
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  # Axes
  axis(side = 1, at = x_axis_seq, cex.axis = 0.8)
  mtext("Time Index", side = 1, line = 2.5, cex = 1)
  axis(side = 2, at = y_axis_seq, cex.axis = 0.5)
  mtext("Diurnally Adjusted Duration", side = 2, line = 2.5, cex = 0.8)
  
  # Add fitted line (black)
  lines(data[[model_name]], col = "black", lty = "solid")
  
  # Model label for legend
  legend_label = switch(model_name,
                        "Wei"  = expression("SCD"["Wei"] * " fitted, " * italic(e)^hat(psi)["i"]),
                        "GG"   = expression("SCD"["GG"] * " fitted, " * italic(e)^hat(psi)["i"]),
                        "Burr" = expression("SCD"["Burr"] * " fitted, " * italic(e)^hat(psi)["i"]),
                        "GB2"  = expression("SCD"["GB2"] * " fitted, " * italic(e)^hat(psi)["i"]),
                        "EGIG" = expression("SCD"["EGIG"] * " fitted, " * italic(e)^hat(psi)["i"]),
                        "TSCD" = expression("TSCD"["EGIG"] * " fitted, " * italic(e)^hat(psi)["i"])
  )
  
  # Legend
  legend("topleft",
         legend = c(expression("Duration, " * italic(x)["i"]),
                    legend_label),
         col = c("grey", "black"), lty = c("solid", "solid"),
         bty = "n", cex = 0.6)
  
  dev.off()
}

# --- 2. Combined plot (grey obs + colored fitted lines for all models) ---
png(filename = "AAPL_obs_and_all_models_fitted.png",
    width = 1400, height = 800, res = 150)

par(mar = c(4, 4, 1, 1), xpd = TRUE)

# Observed durations
plot(data$x.use, type = "l", col = "grey", lty = "solid",
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")

axis(side = 1, at = x_axis_seq)
mtext("Time Index", side = 1, line = 2.5, cex = 1)
axis(side = 2, at = y_axis_seq, cex.axis = 0.8)
mtext("Diurnally Adjusted Duration", side = 2, line = 2.5, cex = 0.8)

# Add fitted lines for all models
for (i in seq_along(models)) {
  lines(data[[models[i]]], col = cols_all[i], lty = ltys_all[i])
}

# Combined legend
legend("topleft",
       legend = c(expression("Observed, " * italic(x)["i"]),
                  expression("SCD"["Wei"] * " fitted, " * italic(e)^hat(psi)["i"]),
                  expression("SCD"["GG"] * " fitted, " * italic(e)^hat(psi)["i"]),
                  expression("SCD"["Burr"] * " fitted, " * italic(e)^hat(psi)["i"]),
                  expression("SCD"["GB2"] * " fitted, " * italic(e)^hat(psi)["i"]),
                  expression("SCD"["EGIG"] * " fitted, " * italic(e)^hat(psi)["i"]),
                  expression("TSCD"["EGIG"] * " fitted, " * italic(e)^hat(psi)["i"])),
       col = c("grey", cols_all),
       lty = c("solid", ltys_all),
       bty = "n", cex = 0.5)

dev.off()
