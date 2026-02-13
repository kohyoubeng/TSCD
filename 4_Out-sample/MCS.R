###############################################################################
# Model Confidence Set (MCS) Forecast Evaluation — SCD & TSCD Models
#
# This script evaluates and compares forecast performance across multiple
# duration models using loss functions and the MCS procedure.
#
# The procedure:
# • Reads observed durations and model forecasts
# • Constructs forecast loss matrices using:
#     – MSFE (Mean Squared Forecast Error)
#     – MAFE (Mean Absolute Forecast Error)
#     – QLIKE loss
# • Saves the loss matrices to CSV files
# • Applies the MCS procedure (Hansen et al.) to identify the superior model set
#
# The SAME script applies to both AAPL and TSLA datasets.
#
# --- To use for TSLA ---
# • Replace input file:
#     "AAPL forc.use, SCDs, TSCD & r.csv"  →  "TSLA forc.use, SCDs, TSCD & r.csv"
# • Replace output filenames beginning with "AAPL" to "TSLA"
#
# Notes:
# • MUST include the following columns in the input file:
#     forc.use, Wei, GG, Burr, GB2, EGIG, TSCD
# • Remaining columns = model forecasts
###############################################################################

rm(list = ls())

# Need to combine manually
# col headings (example): forc.use, Wei, GG, Burr, GB2, EGIG, TSCD, r
data = read.csv("AAPL forc.use, SCDs, TSCD & r.csv", header = TRUE, 
                stringsAsFactors = FALSE)

# ---- sanity checks ----
if (!all(c("forc.use","Wei","GG","Burr","GB2","EGIG","TSCD") %in% colnames(data))) {
  stop("Input file must contain columns: forc.use, Wei, GG, Burr, GB2, EGIG, TSCD")
}

# Prepare containers (not strictly necessary but kept for clarity)
out.MCS_MSFE = NULL
out.MCS_MAFE = NULL
out.MCS_QLIKE = NULL
out.MCS_SSM  = NULL

# True observed values
y = data$forc.use

# Model forecasts only (exclude forc.use and r)
fhat = data[, c("Wei","GG","Burr","GB2","EGIG","TSCD")]

# ---- Loss computations ----
# MSFE
MCS_MSFE = (y - fhat)^2

# MAFE
MCS_MAFE = abs(y - fhat)

# QLIKE: guard against non-positive forecasts or observations
# replace non-positive values with a tiny positive number to avoid -Inf/NaN
eps = .Machine$double.eps * 1e3
fhat_pos = fhat
fhat_pos[fhat_pos <= 0] = eps
y_pos = y
y_pos[y_pos <= 0] = eps

MCS_QLIKE = (y_pos / fhat_pos) - log(y_pos / fhat_pos) - 1

# Save loss matrices to CSV
write.csv(MCS_MSFE,  "AAPL MCS_MSFE.csv",  row.names = FALSE)
write.csv(MCS_MAFE,  "AAPL MCS_MAFE.csv",  row.names = FALSE)
write.csv(MCS_QLIKE, "AAPL MCS_QLIKE.csv", row.names = FALSE)

# ---- Run MCS (drop first row) ----
# We remove the first row ([-1, ]) as discussed so MCS excludes the first/unstable forecast

# MSFE
SSM.MSFE = MCS::MCSprocedure(Loss = MCS_MSFE[-1, ], alpha = 0.05, B = 5000, 
                             statistic = "TR")
sink("AAPL MCS MSFE.txt", append = FALSE)
print(SSM.MSFE)
sink()

# MAFE
SSM.MAFE = MCS::MCSprocedure(Loss = MCS_MAFE[-1, ], alpha = 0.05, B = 5000, 
                             statistic = "TR")
sink("AAPL MCS MAFE.txt", append = FALSE)
print(SSM.MAFE)
sink()

# QLIKE
SSM.QLIKE = MCS::MCSprocedure(Loss = MCS_QLIKE[-1, ], alpha = 0.05, B = 5000, 
                              statistic = "TR")
sink("AAPL MCS QLIKE.txt", append = FALSE)
print(SSM.QLIKE)
sink()
