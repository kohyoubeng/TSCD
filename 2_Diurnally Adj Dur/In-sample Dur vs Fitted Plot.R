###############################################################################
# This script produces the duration and diurnally adjusted duration plots for
# the in-sample period. The plotting procedure is identical for AAPL and TSLA.
# Differences in axis() tick marks and limits reflect different data scales,
# not methodological changes.
#
# --- To reproduce TSLA figures ---
# • Replace input file with: "TSLA Dur Plot (t).csv"
# • Change output filenames accordingly
# • Adjust axis() tick intervals and limits if necessary to reflect TSLA's
#   duration scale. This does not affect the plotting methodology.
#
# All other plotting commands remain unchanged.
###############################################################################

rm(list = ls())
par(mar=c(4,4,1,1), xpd=TRUE)

###############################################################################
# Load In-sample data
data = read.csv("AAPL Dur Plot (t).csv", header=TRUE)

###############################################################################
# === 1) SAVE DURATION PLOT ===
png("AAPL_Duration_plot.png", width = 2000, height = 1200, res = 200)

plot(data$L.dur, type='l', xaxt="n", yaxt="n", xlab="", ylab="")
axis(1, at=seq(0, length(data$L.dur), 500), cex.axis=1.2)
mtext("Time Index", side=1, line=2.5, cex=1.5)

axis(2, at=seq(0, max(data$L.dur), 1000), cex.axis=1.2)
mtext("Observed Duration", side=2, line=2.5, cex=1.5)

dev.off()  # <=== save first PNG

###############################################################################
# === 2) SAVE DIURNALLY ADJUSTED DURATION PLOT ===
png("AAPL_AdjustedDuration_plot.png", width = 2000, height = 1200, res = 200)

plot(data$x.use, type='l', xaxt="n", yaxt="n", xlab="", ylab="")
axis(1, at=seq(0, length(data$L.dur), 500), cex.axis=1.2)
mtext("Time Index", side=1, line=2.5, cex=1.5)

axis(2, at=seq(0, 50, 5), cex.axis=1.2)
mtext("Diurnally adjusted duration", side=2, line=2.5, cex=1.2)

dev.off()  # <=== save second PNG
