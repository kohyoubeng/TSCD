###############################################################################
# Calculate Average Parameters, SD & 95% CI
#
# This script summarizes estimation results from multiple runs:
#   • Computes mean, SD, and 95% CI for model parameters
#   • Computes mean only for fit statistics (LogL, AIC, BIC)
#
# --- To reproduce for TSLA ---
# • Ensure files match pattern "TSLA Est .*SCD.csv"
# • Output summary files will be named like "TSLA Summary <ModelName>.csv"
#
# Notes:
# • MUST get the output CSV files of estimated parameter values for AAPL or TSLA
# • The stock name in file patterns/output is the only difference between AAPL and TSLA runs.
###############################################################################

rm(list = ls())

# List all estimation result files (includes SCD and TSCD models)
files = list.files(pattern = "^AAPL Est .*SCD\\.csv$")

# Function to summarize 500 by LogL
summarize = function(df, n = 500) {
  
  # Keep only numeric columns
  df_num = df[, sapply(df, is.numeric)]
  
  # Identify fit columns and parameter columns
  fit_cols = c("LogL", "AIC", "BIC")
  param_cols = setdiff(colnames(df_num), fit_cols)
  
  # Remove rows with any NA in parameters
  df_clean = df_num[complete.cases(df_num[, param_cols]), ]
  
  # Order by LogL descending and take n rows
  df = df_clean[order(df_clean$LogL, decreasing = TRUE), ][1:n, ]
  
  # Save the 500 clean rows for later inspection
  write.csv(df, paste0("T500_", sub(".csv", "", basename(f)), ".csv"), 
            row.names = FALSE)
  
  # Function to calculate 95% percentile CI
  ci95 = function(x) {
    x = sort(x[!is.na(x)])
    n = length(x)
    if (n == 0) return(c(Lower = NA, Upper = NA))
    lower_idx = ceiling(0.025 * n)
    upper_idx = floor(0.975 * n)
    if (lower_idx < 1) lower_idx = 1
    if (upper_idx < 1) upper_idx = 1
    c(Lower = x[lower_idx], Upper = x[upper_idx])
  }
  
  # Calculate mean, sd, and CI for parameter columns
  mean_vals = sapply(df[, param_cols, drop = FALSE], mean)
  sd_vals = sapply(df[, param_cols, drop = FALSE], sd)
  ci_vals = t(sapply(df[, param_cols, drop = FALSE], ci95))
  
  param_result = data.frame(
    Parameter = names(mean_vals),
    Mean = mean_vals,
    SD = sd_vals,
    CI_Lower = ci_vals[, "Lower"],
    CI_Upper = ci_vals[, "Upper"],
    stringsAsFactors = FALSE
  )
  
  # Calculate mean only for fit statistics (LogL, AIC, BIC)
  fit_mean = colMeans(df[, fit_cols, drop = FALSE])
  fit_result = data.frame(
    Parameter = names(fit_mean),
    Mean = fit_mean,
    SD = NA,
    CI_Lower = NA,
    CI_Upper = NA,
    stringsAsFactors = FALSE
  )
  
  # Combine parameter and fit statistics results
  result = rbind(param_result, fit_result)
  return(result)
}

# Loop through files and save summary results
for (f in files) {
  df = read.csv(f, header = TRUE)
  summary_df = summarize(df, n = 500)
  
  # Generate output file name (replace "Est" with "Summary")
  out_name = sub("Est", "Summary", f)
  
  # Save summary CSV
  write.csv(summary_df, out_name, row.names = FALSE)
  
  cat("✅ Saved summary for:", f, "→", out_name, "\n")
}

cat("🎯 All summaries and T500 subsets saved successfully.\n")
