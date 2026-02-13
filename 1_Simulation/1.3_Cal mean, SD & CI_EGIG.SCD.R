###############################################################################
# This script computes summary statistics for the EGIG.SCD model parameter 
# estimates, including Mean, Bias, Standard Deviation (SD), 
# Mean Squared Error (MSE), and 95% percentile confidence intervals.
#
# The statistics are calculated for different numbers of replications
# (R = 100, 500, 1000) and exported to CSV files for reporting.
###############################################################################

rm(list = ls())

# Read data
# data = read.csv("Sim Est_EGIG.SCD_n=1000 (R=1000).csv", header = TRUE)
# data = read.csv("Sim Est_EGIG.SCD_n=3000 (R=1000).csv", header = TRUE)
data = read.csv("Sim Est_EGIG.SCD_n=6500 (R=1000).csv", header = TRUE)

# True parameter values (must be named to match column names in param_cols)
true.value = c(miu = 0.6000, 
               be  = 0.9700, 
               sig = 0.1000, 
               lam = 0.4000, 
               del = 1.1000, 
               w   = 0.2000)

# Function to calculate 95% CI by sorting (percentile method)
calc_CI = function(x) {
  sorted_x = sort(x)
  n = length(sorted_x)
  lower = sorted_x[ceiling(0.025 * n)]
  upper = sorted_x[floor(0.975 * n)]
  c(lower, upper)
}

# Modified summary stats function
summary_stats = function(df, true_vals) {
  res = lapply(names(df), function(param) {
    x = df[[param]]
    m = mean(x)
    s = sd(x)
    bias = m - true_vals[param]
    mse = mean((x - true_vals[param])^2)
    ci = calc_CI(x)
    
    data.frame(
      Parameter = param,
      Mean = m,
      Bias = bias,
      SD = s,
      MSE = mse,
      CI_lower = ci[1],
      CI_upper = ci[2]
    )
  })
  res_df = do.call(rbind, res)   # proper stacked data frame
  return(res_df)
}

# Parameters only (exclude LogLike)
param_cols = c("miu","be","sig","lam","del","w")

# First 100 rows
res_100 = summary_stats(data[1:100, param_cols], true.value)

# First 500 rows
res_500 = summary_stats(data[1:500, param_cols], true.value)

# First 1000 rows
res_1000 = summary_stats(data[1:1000, param_cols], true.value)

# Combine results into a list
results = list(
  R100 = res_100,
  R500 = res_500,
  R1000 = res_1000
)

# Print results
print(results)

# Save each as CSV
write.csv(res_100, "EGIG.SCD_summary_R100.csv", row.names = FALSE)
write.csv(res_500, "EGIG.SCD_summary_R500.csv", row.names = FALSE)
write.csv(res_1000, "EGIG.SCD_summary_R1000.csv", row.names = FALSE)
