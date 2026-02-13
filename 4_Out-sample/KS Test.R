###############################################################################
# PIT KS Test
#
# This script evaluates the uniformity of PIT values for all SCD and TSCD models
# using the KS test.
#
# The procedure:
# • Performs the KS test for each model against the uniform distribution U(0,1)
# • Captures KS statistics and p-values
# • Saves the overall results to a text file for documentation
#
# The SAME script can be applied to other stocks (e.g., TSLA) with minimal changes.
#
# --- To use for TSLA ---
# • Replace input file:
#     "AAPL All PITs_SCD & TSCD models.csv" → "TSLA All PITs_SCD & TSCD models.csv"
# • Replace output file:
#     "AAPL Overall KS Test for All SCD & TSCD Models.txt" → 
#     "TSLA Overall KS Test for All SCD & TSCD Models.txt"
###############################################################################

rm(list = ls())

#Get x_i forecast and parameter estimates (combine manually)
read.data = read.csv("AAPL All PITs_SCD & TSCD models.csv", header=TRUE)

# Initialize a list to store KS test results for each model
ks_results = list()

# Loop through each column of PIT values (each representing a model)
for (model_name in colnames(read.data)) {
  
  # Extract PIT values for the current model
  PITs = read.data[[model_name]]
  
  # Perform the Kolmogorov-Smirnov test against a uniform distribution
  ks_results[[model_name]] = ks.test(PITs, "punif", 0, 1)
  
}

# Open a sink to capture output in a text file
sink("AAPL Overall KS Test for All SCD & TSCD Models.txt", append = FALSE)

# Display the title
cat("AAPL Overall KS Test Results for All SCD & TSCD Models\n")
cat("=====================================\n\n")

# Display the KS test results
print(ks_results)

# Close the sink
sink()


