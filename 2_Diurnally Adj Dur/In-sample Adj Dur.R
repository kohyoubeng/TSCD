###############################################################################
# This script is used for both AAPL and TSLA datasets. The methodology and all
# processing steps are identical; only the input data file and the event
# threshold (crit) differ.
#
# --- To run for TSLA ---
# 1. Replace the input file with: "TSLA Data_29 Sep-17 Oct 2025.csv"
# 2. Change the event threshold to: crit = 0.19
#
# The remainder of the script is unchanged. Variable names referencing AAPL
# reflect the original dataset used during development and do not affect
# reproducibility.
###############################################################################

rm(list = ls())

# Load raw trade-by-trade AAPL data
AAPL = read.csv("AAPL Data_29 Sep-17 Oct 2025.csv", header=TRUE)

nrow(AAPL)

# ===================== Event Definition =====================
crit = 0.05 #use crit = 0.19 for TSLA

# Initialize event indicator vector
event = numeric(nrow(AAPL))

for(i in 2:nrow(AAPL)){
  # Event if absolute price change ≥ threshold
  if(abs(AAPL$Price[i]-AAPL$Price[i-1])>=crit ){
    
    event[i] = 1
    
  } else {
    
    event[i] = 0
    
  }
}

# Attach event indicator to dataset
AAPLdata = cbind(AAPL, event)

# ===================== Arrival Times =====================
arr.times = AAPLdata$t[event==1]

# Keep first observation + all event observations
L.Cdata = rbind(AAPLdata[1,], AAPLdata[AAPLdata$event==1,])
write.csv(L.Cdata, "Lilian Clean Data.csv")
nrow(L.Cdata)

# ===================== Duration Calculation =====================
L.dur = numeric(nrow(L.Cdata))
for(i in 2:nrow(L.Cdata)){
  L.dur[i] = L.Cdata$t[i] - L.Cdata$t[i-1]
}

# Combine durations with event data
L.com.dur = data.frame(L.Cdata, L.dur)

# Remove zero/negative durations
L.use.dur = rbind(L.com.dur[1,], L.com.dur[L.com.dur$L.dur>0,])
write.csv(L.use.dur, "Lilian All Data with Dur.csv")

# Keep relevant variables
dat.tim.dur = L.use.dur[,c('Date', 't', 'Volume', 'L.dur')]
write.table(dat.tim.dur, file="save date, time, volume, dur.csv", sep=",", row.names = F)

###############################################################################
# ===================== Select Date Range =====================
sel.dat.tim.dur1 = with(dat.tim.dur, dat.tim.dur[(Date < "17-10-2025"), ])

# Convert Date column to proper Date format
dat.tim.dur$Date = as.Date(dat.tim.dur$Date, format = "%d/%m/%Y")

# Subset in-sample window: 29 Sep – 10 Oct 2025
sel.dat.tim.dur1 = subset(dat.tim.dur,
                          Date >= as.Date("2025-09-29") &
                            Date <= as.Date("2025-10-10"))

# Remove first row (no prior duration reference)
sel.dat.tim.dur = sel.dat.tim.dur1[2:nrow(sel.dat.tim.dur1), ]
write.table(sel.dat.tim.dur, file="select date, time, volume & dur.csv", 
            sep=",", row.names = F)

# Extract vectors
t = sel.dat.tim.dur$t
vol = sel.dat.tim.dur$Volume
dur = sel.dat.tim.dur$L.dur
n.sel = length(t)

###############################################################################
# ===================== Diurnal Duration Adjustment =====================
cf_per_day = function(c.date, c.t, c.dur) {
  
  adj.dur = numeric(length(c.t))
  unique_dates = unique(c.date)
  
  for (d in unique_dates) {
    day_index = which(c.date == d)
    t_day  = c.t[day_index]
    dur_day = c.dur[day_index]
    
    # --- Step 1: Winsorize (cap extreme durations) ---
    q99 = quantile(dur_day, 0.99, na.rm = TRUE)
    dur_day[dur_day > q99] = q99
    
    ln.dur = log(dur_day)
    
    # --- Step 2: Construct regressors ---
    f1 = -((t_day - 43200)/14400)^2
    f2 = -((t_day - 48300)/9300)^2
    f3 = ifelse(t_day < 43200, -((t_day - 38700)/7500)^2, 0)
    f4 = ifelse(t_day >= 43200, -((t_day - 48600)/9000)^2, 0)
    f5 = ifelse(t_day >= 34200 & t_day <= 34500, 1, 0)
    f6 = ifelse(t_day > 34500 & t_day <= 34800, 1, 0)
    f7 = ifelse(t_day >= 55800 & t_day <= 57600, 1, 0)
    
    reg.data = data.frame(ln.dur, f1, f2, f3, f4, f5, f6, f7)
    
    # --- Step 3: Robust regression ---
    fit = MASS::rlm(ln.dur ~ f1 + f2 + f3 + f4 + f5 + f6 + f7, data = reg.data)
    
    # --- Step 4: Predict & adjust ---
    d.ti = predict(fit, newdata = reg.data)
    f.ti = exp(d.ti)
    adj_day = dur_day / f.ti
    
    adj.dur[day_index] = adj_day
  }
  
  return(adj.dur)
}


###############################################################################
# Apply per-day duration cleaning
x.use = cf_per_day(sel.dat.tim.dur$Date, t, dur)
out.x.use = cbind(x.use, vol)

write.csv(out.x.use, "Adj dur AAPL In-sample (t)_29 Sep-10 Oct 2025.csv")

###############################################################################
# Final dataset for modeling
final_output = data.frame(
  Date = sel.dat.tim.dur$Date,
  t = t,
  L.dur = dur,
  x.use = x.use
)

write.csv(final_output,
          "AAPL Dur Plot (t).csv", row.names = FALSE)


