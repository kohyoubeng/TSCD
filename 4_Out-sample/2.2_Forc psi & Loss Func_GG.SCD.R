# =============================================================================
# One-Step-Ahead Latent State Forecasting and Loss Func – GG.SCD
#
# NOTE:
# This master code is written using AAPL as the example dataset.
#
# The SAME code structure applies to TSLA without any modification to
# the model or algorithm.
#
# To run for TSLA:
# • Replace "AAPL x.use roll.csv" with "TSLA x.use roll.csv"
# • Replace "AAPL Roll pars_GG.SCD.csv" with "TSLA Roll pars_GG.SCD.csv"
# • Replace all output file names beginning with "AAPL" to "TSLA"
# • Replace "AAPL forc clean obs.csv" with "TSLA forc clean obs.csv"
#
# No changes to the particle filter, forecasting procedure, or
# forecast evaluation metrics are required.
#
# Parameter values used in forecasting are obtained from the rolling 
# estimation output CSV file.
# =============================================================================

rm(list = ls())

# ===================== Load Data =====================
x.read = read.csv("AAPL x.use roll.csv", header=TRUE)
par.read = read.csv("AAPL Roll pars_GG.SCD.csv", header=TRUE)

N = nrow(x.read) #sample size (row)
P = 1000         #number of particles (column)
nf = nrow(par.read)  #number of forecast

# ===================== Allocate memory =====================
psi = matrix(data=0, nrow=N, ncol=P)
psi.est = matrix(data=0, nrow=N, ncol=P)
un.weight = matrix(data=0, nrow=N, ncol=P)
n.weight = matrix(data=0, nrow=N, ncol=P)
resample.location = matrix(data=0, nrow=N, ncol=P)
resamples = matrix(data=0, nrow=N, ncol=P)

fpsi = matrix(data=0, nrow=nf, ncol=P)  # store forecast psi at each step

forc1 = numeric(nf)

# ===================== Particle Filter Function =====================
my_func = function(x.obs, miu, be, sig, a, p){
  
  for(i in 1:N){
    if(i == 1){
      
      # -------- Initial latent state from stationary distribution --------
      psi0 = rnorm(n=P, mean=miu, sd = sig/sqrt(1-be^2))
      
      #Generate state from model:
      psi[i, ] = miu*(1-be) + be*psi0 + rnorm(P, mean=0, sd=sig)
      
    } else {
      
      #Generate state from model:
      psi[i, ] = miu*(1-be) + be*psi.est[i-1, ] + rnorm(P, mean=0, sd=sig)
      
    }
    
    # ================= Calculate weight =================
    #Compute unnormalized weight:
    g1 = gamma(p)
    g2 = gamma(p + 1/a)
    
    un.weight[i, ] = a * g2^(a*p) / g1^(a*p + 1) *
      x.obs[i]^(a*p - 1) / (exp(psi[i, ]))^(a*p) *
      exp(-((x.obs[i] * g2) / (exp(psi[i, ]) * g1))^a)
    
    #Compute normalized weight:
    n.weight[i, ] = un.weight[i, ] / sum(un.weight[i, ])
    
    # ================= Resampling (SIR) =================
    resample.location[i, ] = sample(1:P, size = P, replace = TRUE, 
                                    prob = n.weight[i, ])
    psi.est[i, ] = psi[i, resample.location[i, ]]
    
  }
  
  return(list(psi=psi, psi.est=psi.est))
  
}

# ===================== Rolling Forecast Loop =====================
for(j in 1:nf) {
  
  print(j)
  
  #Read data
  x.obs = x.read[, j]
  
  # Load estimated parameters for this window
  miu = par.read$miu[j]
  be = par.read$be[j]
  sig = par.read$sig[j]
  a = par.read$a[j]
  p = par.read$p[j]
  
  # Run particle filter on this window
  res = my_func(x.obs, miu, be, sig, a, p)
  
  # Last filtered state distribution ψ_n
  psi.prev = res$psi.est[N, ]
  
  # ---- FORECAST LATENT STATE ψ_{n+1} ----
  fpsi[j, ] = miu*(1-be) + be*psi.prev + rnorm(P, mean=0, sd=sig)
  
  # ---- FORECAST OBSERVATION ----
  forc1[j] = mean(exp(fpsi[j, ]))  # E[exp(psi_{n+1})]
}

# Save forecast latent states
colnames(fpsi) = paste0("Parti", 1:P)
rownames(fpsi) = paste0("fpsi", 1:nf)
write.csv(fpsi, "AAPL fpsi_GG.SCD.csv", row.names=TRUE)
write.csv(forc1, "AAPL forc_GG.SCD.csv", row.names=FALSE)

################################################################################
# ===================== Forecast Evaluation =====================
################################################################################

# Load actual observed future durations
fobs.read = read.csv("AAPL forc clean obs.csv", header=TRUE)
fobs = array(unlist(fobs.read$forc.use))

# Reload forecasts
forecast.df = read.csv("AAPL forc_GG.SCD.csv", header=TRUE)

forc1 = forecast.df$x

# ---- Forecast Accuracy Metrics ----
MSFE_1 = mean((fobs - forc1)^2)
MAFE_1 = mean(abs(fobs - forc1))
QLIKE_1 = mean((fobs / forc1) - log(fobs / forc1) - 1, na.rm = TRUE)

Out_1 = rbind(MSFE_1, MAFE_1, QLIKE_1)
colnames(Out_1) = "AAPL GG.SCD"

# Combine into one data frame
Out_combined = data.frame(
  Metric = c("MSFE", "MAFE", "QLIKE"),
  forc1 = Out_1[, 1]
)

# Save to CSV
write.csv(Out_combined, "AAPL Forc Error_GG.SCD.csv", row.names = FALSE)

Out_1



