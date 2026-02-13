###############################################################################
# DM Test
#
# This script compares forecast accuracy of multiple models for AAPL using
# the DM Test.
#
# The procedure:
# • Reads observed durations and model forecasts from CSV
# • Computes forecast errors using three loss functions:
#     – MSFE  (Mean Squared Forecast Error)
#     – MAFE  (Mean Absolute Forecast Error)
#     – QLIKE (Quasi-Likelihood Loss)
# • Calculates pairwise loss differences for all models
# • Computes long-run variance and DM statistics for each model pair
# • Calculates p-values and determines statistical significance
# • Stores results in both matrix and tabular CSV files
#
# The SAME script can be applied to other stocks (e.g., TSLA) with minimal changes.
#
# --- To use for TSLA ---
# • Replace input file:
#     "AAPL forc.use, SCDs, TSCD & r.csv" → "TSLA forc.use, SCDs, TSCD & r.csv"
# • Replace output filenames beginning with "AAPL" to "TSLA"
###############################################################################

rm(list = ls())

# Load needed library
library(sandwich)

# Read data
data = read.csv("AAPL forc.use, SCDs, TSCD & r.csv", header = TRUE)

x = data$forc.use             # actual values
x_hat = data[, 2:7]  # forecasts for all models

h = 1  # forecast horizon
num_models = ncol(x_hat)
model_names = colnames(x_hat)

# Loss functions
msfe_loss = function(y, f) {
  return((y - f)^2)
}

mafe_loss = function(y, f) {
  return(abs(y - f))
}

qlike_loss = function(y, f) {
  y_div_f = y / f
  return(y_div_f - log(y_div_f) - 1)
}

############################################################################
# Custom function to compute long-run variance (Diebold-Mariano style)
dm_variance = function(d, h = 1) {
  T = length(d)
  d_mean = mean(d)
  centered_d = d - d_mean
  gamma_0 = sum(centered_d^2) / T
  
  if (h == 1) {
    return(gamma_0)
  } else {
    gamma_sum = 0
    for (k in 1:(h - 1)) {
      weight = 1 - (k / h)
      gamma_k = sum(centered_d[(k + 1):T] * centered_d[1:(T - k)]) / T
      gamma_sum = gamma_sum + 2 * weight * gamma_k
    }
    return(gamma_0 + gamma_sum)
  }
}

# DM test using DM (1995) long-run variance###############################
dm_test_custom = function(loss_diff, h = 1) {
  d_bar = mean(loss_diff)
  T = length(loss_diff)
  var_d = dm_variance(loss_diff, h)
  DM_stat = d_bar / sqrt(var_d / T)
  p_val = 2 * (1 - pnorm(abs(DM_stat)))
  return(list(statistic = DM_stat, p.value = p_val))
}

# Initialize matrices to store p-values
pval_MSFE = matrix(NA, nrow = num_models, ncol = num_models)
pval_MAFE = matrix(NA, nrow = num_models, ncol = num_models)
pval_QLIKE = matrix(NA, nrow = num_models, ncol = num_models)

rownames(pval_MSFE) = colnames(pval_MSFE) = model_names
rownames(pval_MAFE) = colnames(pval_MAFE) = model_names
rownames(pval_QLIKE) = colnames(pval_QLIKE) = model_names

# Data frame to store results
dm_df = data.frame(Model_1 = character(),
                   Model_2 = character(),
                   Loss_Function = character(),
                   P_Value = numeric(),
                   stringsAsFactors = FALSE)

# Loop over model pairs
for (i in 1:(num_models - 1)) {
  for (j in (i + 1):num_models) {
    
    d_msfe = msfe_loss(x, x_hat[, i]) - msfe_loss(x, x_hat[, j])
    d_mafe = mafe_loss(x, x_hat[, i]) - mafe_loss(x, x_hat[, j])
    d_qlike = qlike_loss(x, x_hat[, i]) - qlike_loss(x, x_hat[, j])
    
    dm_msfe = dm_test_custom(d_msfe, h = h)
    dm_mafe = dm_test_custom(d_mafe, h = h)
    dm_qlike = dm_test_custom(d_qlike, h = h)
    
    dm_df = rbind(dm_df, data.frame(Model_1 = model_names[i],
                                    Model_2 = model_names[j],
                                    Loss_Function = "MSFE",
                                    P_Value = round(dm_msfe$p.value, 4),
                                    stringsAsFactors = FALSE))
    
    dm_df = rbind(dm_df, data.frame(Model_1 = model_names[i],
                                    Model_2 = model_names[j],
                                    Loss_Function = "MAFE",
                                    P_Value = round(dm_mafe$p.value, 4),
                                    stringsAsFactors = FALSE))
    
    dm_df = rbind(dm_df, data.frame(Model_1 = model_names[i],
                                    Model_2 = model_names[j],
                                    Loss_Function = "QLIKE",
                                    P_Value = round(dm_qlike$p.value, 4),
                                    stringsAsFactors = FALSE))
    
    pval_MSFE[i, j] = round(dm_msfe$p.value, 4)
    pval_MAFE[i, j] = round(dm_mafe$p.value, 4)
    pval_QLIKE[i, j] = round(dm_qlike$p.value, 4)
  }
}

dm_df$Significance = ifelse(dm_df$P_Value < 0.05,
                            "Significant difference",
                            "Not significant difference")

write.csv(pval_MSFE, "AAPL DM MSFE matrix_Paper2.csv")
write.csv(pval_MAFE, "AAPL DM MAFE matrix_Paper2.csv")
write.csv(pval_QLIKE, "AAPL DM QLIKE matrix_Paper2.csv")

write.csv(dm_df, "AAPL DM results_Paper2.csv", row.names = FALSE)

print(dm_df)

