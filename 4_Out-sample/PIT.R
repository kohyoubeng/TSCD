###############################################################################
# PIT
#
# This script produces the analyses of PIT values for all SCD and TSCD models.
#
# The procedure:
# • Get PIT values
# • Generates three types of plots for each model:
#     1. Scatter plot of PIT values vs. time index
#     2. Histogram of PIT values with 95% confidence intervals for uniformity
#     3. Empirical CDF vs. theoretical CDF (U(0,1)) with a 45-degree reference line
# • Saves all plots as PNG files with descriptive filenames
# • Outputs status messages to track progress for each model
#
# The SAME script can be applied to other stocks (e.g., TSLA) with minimal changes.
#
# --- To use for TSLA ---
# • Replace input file:
#     "AAPL All PITs_SCD & TSCD models.csv" → "TSLA All PITs_SCD & TSCD models.csv"
# • Replace output filenames in `png()` calls:
#     "AAPL PIT Scatter_", "AAPL PIT Histogram_", "AAPL PIT CDF_" → 
#     "TSLA PIT Scatter_", "TSLA PIT Histogram_", "TSLA PIT CDF_"
###############################################################################

rm(list = ls())

read.data = read.csv("AAPL All PITs_SCD & TSCD models.csv", header = TRUE)

# Get all PIT columns starting with "PIT."
models_full = grep("^PIT\\.", colnames(read.data), value = TRUE)

m = 10  # histogram bins

for (model_col in models_full) {
  u = read.data[[model_col]]
  u = na.omit(u)
  
  # Clip PIT values to [0,1] just in case
  u = pmax(pmin(u, 1), 0)
  
  # Extract file name from column header (remove "PIT.")
  model_file_name = sub("^PIT\\.", "", model_col)  # e.g. "Wei.SCD"
  
  # Extract suffix (SCD or TSCD) and distribution name
  # Assuming suffix is last after dot, e.g. Wei.SCD or EGIG.TSCD
  parts = unlist(strsplit(model_file_name, "\\."))
  dist = paste(parts[-length(parts)], collapse = ".")  # distribution name, e.g. "Wei" or "EGIG"
  suffix = parts[length(parts)]                        # "SCD" or "TSCD"
  
  # Scatter plot
  png(filename = paste0("AAPL PIT Scatter_", model_file_name, ".png"), 
      width = 1600, height = 1000, res = 200)
  par(mar = c(5, 5, 4, 2) + 0.1)  # increase left margin to avoid cutting y-axis label
  plot(u, ylab = "PIT Values", xlab = "Time Index", pch = 16, col = "blue", 
       main = "",
       cex.lab = 2,    # make axis labels 2 times default size
       cex.axis = 1.5)   # make tick labels 1.5 times default size
  dev.off()
  
  # Histogram
  n = length(u)
  expected_count = n / m
  se_bin = sqrt(expected_count * (1 - 1/m))
  ci_lower = expected_count - 1.96 * se_bin
  ci_upper = expected_count + 1.96 * se_bin
  
  breaks_vec = seq(0, 1, length.out = m + 1)
  
  png(filename = paste0("AAPL PIT Histogram_", model_file_name, ".png"), 
      width = 1600, height = 1000, res = 200)
  par(mar = c(5, 5, 4, 2) + 0.1)  # increase left margin here too
  hist(u, breaks = breaks_vec, col = "lightblue", border = "darkblue",
       xlab = "PIT Values", ylab = "Frequency", xlim = c(0, 1), main = "",
       cex.lab = 2,    # bigger axis labels
       cex.axis = 1.5)   # bigger tick labels
  abline(h = ci_lower, col = "black", lty = 2, lwd = 2)
  abline(h = ci_upper, col = "black", lty = 2, lwd = 2)
  dev.off()
  
  # CDF plot
  sorted_u = sort(u)
  n.f = length(u)
  empirical_cdf = (1:n.f) / n.f
  theoretical_cdf = sorted_u
  
  png(filename = paste0("AAPL PIT CDF_", model_file_name, ".png"),
      width = 1600, height = 1000, res = 200)
  
  par(mar = c(4.5, 5, 2.5, 2))
  
  plot(theoretical_cdf, empirical_cdf, type = "s",
       col = "red", lwd = 2.8,        # <<< thicker solid line
       xlab = "PIT Values", ylab = "Cumulative Probability",
       xaxp = c(0, 1, 10), yaxp = c(0, 1, 10),
       main = "",
       cex.lab = 1.6,     # axis label size
       cex.axis = 1.3,    # tick label size
       cex.main = 1.8)
  
  # 45-degree theoretical line
  abline(0, 1, col = "blue", lty = 2, lwd = 2.8)   # <<< thicker dotted line
  
  legend(
    x = 0.01, y = 0.98,
    legend = c(
      bquote("Empirical CDF, " * .(suffix)[.(dist)]),
      "Theoretical CDF, U(0,1)"
    ),
    col = c("red", "blue"),
    lty = c(1, 2),
    lwd = c(2.8, 2.8),        # <<< thicker legend lines too
    bty = "o",
    cex = 1.2                   # <<< Larger legend wording
  )
  
  dev.off()
  
  cat("✅ Finished plots for model:", model_file_name, "\n")
}
