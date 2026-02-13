###############################################################################
# Run Estimation Scripts Multiple Times
#
# This script repeats a selected estimation procedure multiple times (e.g., 500 runs):
#   • Loops over run numbers
#   • Sources of estimation script and data **MUST** save in the correct working directory
#   • Optionally removes the previous iteration's results
#   • Saves each run's output with a unique filename including the run number
#
# --- To reproduce for TSLA ---
# • Output files will be named like "TSLA Est_<ModelName>_<RunNumber>.csv"
#
# Notes:
# • The loop count, source script, and file naming are the only differences
#   between AAPL and TSLA runs. All other code is identical.
###############################################################################

for (run in 1:500) {   # repeat 500 times
  cat("Run number:", run, "\n")
  
  # --- optionally remove only results from previous iteration ---
  if (exists("results")) rm(results)
  
  # run the estimation script
  # source("1.1_Est_Wei.SCD.R")
  # source("2.1_Est_GG.SCD.R")
  # source("3.1_Est_Burr.SCD.R")
  # source("4.1_Est_GB2.SCD.R")
  # source("5.1_Est_EGIG.SCD.R")
  source("6.1_Est_EGIG.TSCD.R")
  
  # --- save results uniquely ---
  # out_file = paste0("AAPL Est_Wei.SCD_", run, ".csv")
  # out_file = paste0("AAPL Est_GG.SCD_", run, ".csv")
  # out_file = paste0("AAPL Est_Burr.SCD_", run, ".csv")
  # out_file = paste0("AAPL Est_GB2.SCD_", run, ".csv")
  # out_file = paste0("AAPL Est_EGIG.SCD_", run, ".csv")
  out_file = paste0("AAPL Est_EGIG.TSCD_", run, ".csv")
  
  write.table(results, file = out_file, sep = ",", row.names = FALSE)
}
