###############################################################################
# This script is used for both AAPL and TSLA datasets. The methodology and all
# processing steps are identical; only the input data file and the event
# threshold (crit) differ.
#
# --- To run for TSLA ---
# 1. Replace the input file with: "TSLA 29 Sep-17 Oct (date, t, price, volume).csv"
# 2. Change the event threshold to: crit = 0.19
#
# The remainder of the script is unchanged. Variable names referencing AAPL
# reflect the original dataset used during development and do not affect
# reproducibility.
###############################################################################

rm(list = ls())
library(MASS)

# --------------------------- 1. Read data ------------------------------------
AAPL = read.csv("AAPL 29 Sep-17 Oct (date, t, price, volume).csv", header = TRUE)

# --------------------------- 2. Define event variable -------------------------
crit = 0.05  #use crit = 0.19 for TSLA

event = numeric(nrow(AAPL))

for (i in 2:nrow(AAPL)) {
  if (abs(AAPL$Price[i] - AAPL$Price[i - 1]) >= crit) {
    event[i] = 1
  } else {
    event[i] = 0
  }
}

AAPLdata = cbind(AAPL, event)

# --------------------------- 3. Select event rows -----------------------------
arr.times = AAPLdata$t[event == 1]
L.Cdata = rbind(AAPLdata[1, ], AAPLdata[AAPLdata$event == 1, ])
write.csv(L.Cdata, "Lilian Clean Data.csv", row.names = FALSE)

# --------------------------- 4. Compute durations -----------------------------
L.dur = numeric(nrow(L.Cdata))
for (i in 2:nrow(L.Cdata)) {
  L.dur[i] = L.Cdata$t[i] - L.Cdata$t[i - 1]
}

L.com.dur = data.frame(L.Cdata, L.dur)
L.use.dur = rbind(L.com.dur[1, ], L.com.dur[L.com.dur$L.dur > 0, ])
write.csv(L.use.dur, "Lilian All Data with Dur.csv", row.names = FALSE)

dat.tim.dur = L.use.dur[, c('Date', 't', 'L.dur')]
write.table(dat.tim.dur, file = "save date, time, dur.csv", sep = ",", row.names = FALSE)

# --------------------------- 5. Convert Date ----------------------------------
dat.tim.dur$Date = as.Date(dat.tim.dur$Date, format = "%d/%m/%Y")

# --------------------------- 6. Select samples --------------------------------
sel.dat.tim.dur.in = subset(dat.tim.dur, Date >= as.Date("2025-09-29") & 
                              Date <= as.Date("2025-10-10"))
sel.dat.tim.dur.in = sel.dat.tim.dur.in[2:nrow(sel.dat.tim.dur.in), ]#drop first

sel.dat.tim.dur.out = subset(dat.tim.dur, Date >= as.Date("2025-10-13") & 
                               Date <= as.Date("2025-10-17"))

write.table(sel.dat.tim.dur.in, file = "In-sample.csv", sep = ",", 
            row.names = FALSE)
write.table(sel.dat.tim.dur.out, file = "Out-sample.csv", sep = ",", 
            row.names = FALSE)

###############################################################################
# Function: per-day diurnal adjustment Φ(ti)
###############################################################################

cf_per_day = function(c.date, c.t, c.dur) {
  adj.dur = numeric(length(c.t))
  unique_dates = unique(c.date)
  
  for (d in unique_dates) {
    day_index = which(c.date == d)
    t_day = c.t[day_index]
    dur_day = c.dur[day_index]
    
    # --- Step 1: Winsorize ---
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
    
    # --- Step 3: Robust regression with fallback ---
    fit = tryCatch({
      MASS::rlm(ln.dur ~ f1 + f2 + f3 + f4 + f5 + f6 + f7, data = reg.data)
    }, error = function(e) {
      lm(ln.dur ~ f1 + f2 + f3 + f4 + f5 + f6 + f7, data = reg.data)
    })
    
    # --- Step 4: Predict & adjust ---
    d.ti = predict(fit, newdata = reg.data)
    f.ti = exp(d.ti)
    adj_day = dur_day / f.ti
    
    adj.dur[day_index] = adj_day
  }
  return(adj.dur)
}

###############################################################################
# Apply adjustment to in-sample
###############################################################################

t = sel.dat.tim.dur.in$t
dur = sel.dat.tim.dur.in$L.dur
x.use = cf_per_day(sel.dat.tim.dur.in$Date, t, dur)
write.csv(x.use, "AAPL Adj dur In-sample.csv", row.names = FALSE)

###############################################################################
# Rolling-window construction
###############################################################################

n.in = nrow(sel.dat.tim.dur.in)
n.out = nrow(sel.dat.tim.dur.out)

out.t.rol = list()
out.dur.rol = list()
out.x.use.rol = list()

for (k in 1:n.out) {
  # Rolling window: drop first k in-sample rows, add first k out-sample rows
  sel.dat.tim.dur.rol = rbind(sel.dat.tim.dur.in[-(1:k), ], 
                              sel.dat.tim.dur.out[1:k, ])
  
  # Ensure exact same length as in-sample (n.in)
  if (nrow(sel.dat.tim.dur.rol) > n.in) {
    sel.dat.tim.dur.rol = sel.dat.tim.dur.rol[1:n.in, ]
  }
  
  t.rol = sel.dat.tim.dur.rol$t
  dur.rol = sel.dat.tim.dur.rol$L.dur
  
  out.t.rol[[k]] = t.rol
  out.dur.rol[[k]] = dur.rol
  out.x.use.rol[[k]] = cf_per_day(sel.dat.tim.dur.rol$Date, t.rol, dur.rol)
}

###############################################################################
# Safe combine
###############################################################################

pad_to_len = function(x, len) { length(x) <- len; x }

maxlen = n.in
x.use.padded = pad_to_len(x.use, maxlen)

out.t.rol.padded = lapply(out.t.rol, pad_to_len, len = maxlen)
out.dur.rol.padded = lapply(out.dur.rol, pad_to_len, len = maxlen)
out.x.use.rol.padded = lapply(out.x.use.rol, pad_to_len, len = maxlen)

out.t.final = cbind(x.use.padded, do.call(cbind, out.t.rol.padded))
out.dur.final = cbind(dur, do.call(cbind, out.dur.rol.padded))
out.x.final = cbind(x.use.padded, do.call(cbind, out.x.use.rol.padded))

# ---------- Check expected vs actual ----------
expected_cols = 1 + length(out.x.use.rol)   # 1 baseline + number of roll windows
cat("Expected columns (1 + length(out.x.use.rol)):", expected_cols, "\n")
cat("Actual columns (out.x.final):", ncol(out.x.final), "\n")

# ---------- Name columns sequentially and save ----------
colnames(out.t.final)  = c("t.in", paste0("t.roll", 1:(ncol(out.t.final)-1)))
colnames(out.dur.final) = c("Dur.in", paste0("Dur.roll", 1:(ncol(out.dur.final)-1)))
colnames(out.x.final)  = c("x.in", paste0("x.roll", 1:(ncol(out.x.final)-1)))

write.csv(out.t.final, "AAPL t roll.csv", row.names = FALSE)
write.csv(out.dur.final, "AAPL Dur roll.csv", row.names = FALSE)
write.csv(out.x.final, "AAPL x.use roll.csv", row.names = FALSE)

