###############################################################################
# This script computes summary statistics for the EGIG.TSCD model parameter 
# estimates, including Mean, Bias, Standard Deviation (SD), 
# Mean Squared Error (MSE), and 95% percentile confidence intervals.
#
# The statistics are calculated for different numbers of replications
# (R = 100, 500, 1000) and exported to CSV files for reporting.
###############################################################################

rm(list = ls())

# Read data
data = read.csv("Sim Est_EGIG.TSCD_n=1000 (R=1000).csv", header = TRUE)
# data = read.csv("Sim Est_EGIG.TSCD_n=3000 (R=1000).csv", header = TRUE)
# data = read.csv("Sim Est_EGIG.TSCD_n=6500 (R=1000).csv", header = TRUE)

# True parameter values (must match column names in param_cols)
true.value = c(
  miu1 = 0.3,
  miu2 = 1.2,
  be1 = 0.9,
  be2 = 0.95,
  lam1 = 0.3,
  lam2 = 0.5,
  del1 = 1.05,
  del2 = 1.3,
  w1 = 0.3,
  w2 = 0.1,
  sig1 = 0.15,
  sig2 = 0.06
)

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
  stopifnot(all(names(df) %in% names(true_vals)))  # safety check
  
  res = lapply(names(df), function(param) {
    x = df[[param]]
    m = mean(x)
    s = sd(x)
    bias = m - true_vals[param]
    mse = mean((x - true_vals[param])^2)
    ci = calc_CI(x)
    
    data.frame(
      Parameter = param,
      True_Value = true_vals[param],
      Mean = m,
      Bias = bias,
      SD = s,
      MSE = mse,
      CI_lower = ci[1],
      CI_upper = ci[2]
    )
  })
  
  res_df = do.call(rbind, res)
  return(res_df)
}

# Correct parameter list
param_cols = c("miu1","miu2","be1","be2","sig1","sig2","lam1","lam2",
               "del1","del2","w1","w2")

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
write.csv(res_100, "EGIG.TSCD_summary_R100.csv", row.names = FALSE)
write.csv(res_500, "EGIG.TSCD_summary_R500.csv", row.names = FALSE)
write.csv(res_1000, "EGIG.TSCD_summary_R1000.csv", row.names = FALSE)
