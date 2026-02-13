###############################################################################
# Combine All Valid Output Files
#
# This script loops through all estimation output files for a given
# asset and combines the valid results into a single CSV file. It performs:
#   • Skipping missing files
#   • Skipping files containing NA
#   • Adding a 'Run' column to indicate the simulation/estimation number
#   • Combining all valid files into one data frame
#   • Saving the combined results to a CSV
#
# --- To reproduce TSLA results ---
# • Set pattern = "TSLA Est_<ModelName>_"
# • Output file example: "TSLA Est Wei.SCD combined.csv"
#
# All other steps remain identical between assets. Differences only arise
# from the input filenames corresponding to the asset and model.
###############################################################################

library(dplyr)

# Define pattern (adjust if you use other model names)
# pattern = "AAPL Est_Wei.SCD_"
# pattern = "AAPL Est_GG.SCD_"
# pattern = "AAPL Est_Burr.SCD_"
# pattern = "AAPL Est_GB2.SCD_"
# pattern = "AAPL Est_EGIG.SCD_"
pattern = "AAPL Est_EGIG.TSCD_"

# Initialize an empty list to store valid data
combined_list = list()

# Loop through expected run numbers
for (run in 1:200) {
  file_name = paste0(pattern, run, ".csv")
  
  # Skip if file does not exist
  if (!file.exists(file_name)) {
    cat("Skipping missing file:", file_name, "\n")
    next
  }
  
  # Read file
  df = tryCatch(read.csv(file_name, header = TRUE), error = function(e) NULL)
  
  # Skip if read failed or if it contains NA
  if (is.null(df)) {
    cat("Error reading:", file_name, "\n")
    next
  }
  
  if (any(is.na(df))) {
    cat("Skipping file with NA:", file_name, "\n")
    next
  }
  
  # Add a run identifier column
  df$Run = run
  
  # Store valid data
  combined_list[[length(combined_list) + 1]] = df
}

# Combine all into one data frame
if (length(combined_list) > 0) {
  combined_results = bind_rows(combined_list)
  
  # Save combined results
  # write.csv(combined_results, "AAPL Est Wei.SCD combined.csv", 
  #           row.names = FALSE)
  # write.csv(combined_results, "AAPL Est GG.SCD combined.csv", 
  #           row.names = FALSE)
  # write.csv(combined_results, "AAPL Est Burr.SCD combined.csv", 
  #           row.names = FALSE)
  # write.csv(combined_results, "AAPL Est GB2.SCD combined.csv", 
  #           row.names = FALSE)
  # write.csv(combined_results, "AAPL Est EGIG.SCD combined.csv", 
  #           row.names = FALSE)
  write.csv(combined_results, "AAPL Est EGIG.TSCD combined.csv",
            row.names = FALSE)
  
  # cat("✅ Combined file saved as 'AAPL Est Wei.SCD combined.csv'\n")  
  # cat("✅ Combined file saved as 'AAPL Est GG.SCD combined.csv'\n")  
  # cat("✅ Combined file saved as 'AAPL Est Burr.SCD combined.csv'\n")  
  # cat("✅ Combined file saved as 'AAPL Est GB2.SCD combined.csv'\n")  
  # cat("✅ Combined file saved as 'AAPL Est EGIG.SCD combined.csv'\n")  
  cat("✅ Combined file saved as 'AAPL Est EGIG.TSCD combined.csv'\n")
} else {
  cat("⚠️ No valid files found.\n")
}



