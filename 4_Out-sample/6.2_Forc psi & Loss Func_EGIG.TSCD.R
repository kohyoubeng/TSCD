# =============================================================================
# One-Step-Ahead Latent State Forecasting and Loss Func – EGIG.TSCD
#
# NOTE:
# This master code is written using AAPL as the example dataset.
#
# The SAME code structure applies to TSLA without any modification to
# the model or algorithm.
#
# To run for TSLA:
# • Replace "AAPL x.use roll.csv" with "TSLA x.use roll.csv"
# • Replace "AAPL Roll pars_EGIG.TSCD.csv" with "TSLA Roll pars_EGIG.TSCD.csv"
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
par.read = read.csv("AAPL Roll pars_EGIG.TSCD.csv", header=TRUE)

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

r = numeric(nf)
forc1 = numeric(nf)

# ===================== Particle Filter Function =====================
my_func = function(x.obs, r, miu1, miu2, be1, be2, sig1, sig2, lam1, lam2,
                   del1, del2, w1, w2){
  
  for(i in 1:N){
    if(i == 1){
      if(x.obs[i] <= r){
        
        # -------- Initial latent state from stationary distribution --------
        psi0 = rnorm(n=P, mean=miu1, sd=sig1/sqrt(1-be1^2))
        
        #Generate state from model:
        psi[i, ] = miu1*(1-be1) + be1*psi0 + rnorm(P, mean=0, sd=sig1)
        
        # ================= Calculate weight =================
        #Compute normalized weight:
        b1 = besselK(w1, lam1/del1, expon.scaled = FALSE)
        b2 = besselK(w1, (lam1+1)/del1, expon.scaled = FALSE)
        
        un.weight[i, ] = del1/2 * (b2^lam1 / b1^(lam1+1)) *
          (x.obs[i]^(lam1-1) / exp(psi[i, ])^lam1) *
          exp(-w1/2 * (((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^del1 +
                         ((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^-del1))
        
      } else {
        
        #Give samples:
        psi0 = rnorm(n=P, mean=miu2, sd=sig2/sqrt(1-be2^2))
        
        #Generate state from model:
        psi[i, ] = miu2*(1-be2) + be2*psi0 + rnorm(P, mean=0, sd=sig2)
        
        #Compute normalized weight:
        b1 = besselK(w2, lam2/del2, expon.scaled = FALSE)
        b2 = besselK(w2, (lam2+1)/del2, expon.scaled = FALSE)
        
        un.weight[i, ] = del2/2 * (b2^lam2 / b1^(lam2+1)) *
          (x.obs[i]^(lam2-1) / exp(psi[i, ])^lam2) *
          exp(-w2/2 * (((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^del2 +
                         ((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^-del2))
        
      }
      
    } else {
      
      if(x.obs[i-1] <= r){
        
        #Generate state from model:
        psi[i, ] = miu1*(1-be1) + be1*psi.est[i-1, ] + rnorm(P, mean=0, sd=sig1)
        
        #Compute normalized weight:
        b1 = besselK(w1, lam1/del1, expon.scaled = FALSE)
        b2 = besselK(w1, (lam1+1)/del1, expon.scaled = FALSE)
        
        un.weight[i, ] = del1/2 * (b2^lam1 / b1^(lam1+1)) *
          (x.obs[i]^(lam1-1) / exp(psi[i, ])^lam1) *
          exp(-w1/2 * (((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^del1 +
                         ((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^-del1))
        
      } else {
        
        #Generate state from model:
        psi[i, ] = miu2*(1-be2) + be2*psi.est[i-1, ] + rnorm(P, mean=0, sd=sig2)
        
        #Compute normalized weight:
        b1 = besselK(w2, lam2/del2, expon.scaled = FALSE)
        b2 = besselK(w2, (lam2+1)/del2, expon.scaled = FALSE)
        
        un.weight[i, ] = del2/2 * (b2^lam2 / b1^(lam2+1)) *
          (x.obs[i]^(lam2-1) / exp(psi[i, ])^lam2) *
          exp(-w2/2 * (((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^del2 +
                         ((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^-del2))
        
      }
    }
    
    #Compute normalized weight
    # Remove impossible values
    un.weight[i, ][!is.finite(un.weight[i, ])] = 0
    
    sumw = sum(un.weight[i, ])
    
    # If all weights collapsed → use equal weights
    if(sumw <= 0 || is.na(sumw)){
      n.weight[i, ] = rep(1/P, P)
    } else {
      n.weight[i, ] = un.weight[i, ] / sumw
    }
    
    
    # ================= Resampling (SIR) =================
    resample.location[i, ] = sample(1:P, size = P, replace = TRUE, 
                                    prob = n.weight[i, ])
    psi.est[i, ] = psi[i, resample.location[i, ]]
  }
  
  return(list(psi=psi, psi.est=psi.est))
}

# ===================== Rolling Forecast Loop =====================
for(j in 1:nf){
  
  print(j)
  
  #Read data
  x.obs = x.read[, j]
  
  #threshold value
  r[j] = median(x.obs)
  
  # Load estimated parameters for this window
  miu1 = par.read$miu1[j]
  miu2 = par.read$miu2[j]
  be1 = par.read$be1[j]
  be2 = par.read$be2[j]
  sig1 = par.read$sig1[j]
  sig2 = par.read$sig2[j]
  lam1 = par.read$lam1[j]
  lam2 = par.read$lam2[j]
  del1 = par.read$del1[j]
  del2 = par.read$del2[j]
  w1 = par.read$w1[j]
  w2 = par.read$w2[j]
  
  # Run particle filter on this window
  res = my_func(x.obs, r[j], miu1, miu2, be1, be2, sig1, sig2, lam1, lam2,
                del1, del2, w1, w2)
  
  # Last filtered state distribution ψ_n
  psi.prev = res$psi.est[N, ]
  
  # ---- FORECAST LATENT STATE ψ_{n+1} ----
  if(x.obs[N] <= r[j]){
    
    #forecast psi
    fpsi[j, ] = miu1*(1-be1) + be1*psi.prev + rnorm(P, mean=0, sd=sig1)
    
  } else {
    
    fpsi[j, ] = miu2*(1-be2) + be2*psi.prev + rnorm(P, mean=0, sd=sig2)
    
  }
  
  # ---- FORECAST OBSERVATION ----
  forc1[j] = mean(exp(fpsi[j, ]))  # E[exp(psi_{n+1})]
  
}

# Save forecast latent states
colnames(fpsi) = paste0("Parti", 1:P)
rownames(fpsi) = paste0("fpsi", 1:nf)
write.csv(fpsi, "AAPL fpsi_EGIG.TSCD.csv", row.names=TRUE)
write.csv(r, "AAPL threshold r_EGIG.TSCD.csv", row.names=FALSE)
write.csv(forc1, "AAPL forc_EGIG.TSCD.csv", row.names=FALSE)

################################################################################
# ===================== Forecast Evaluation =====================
################################################################################

# Load actual observed future durations
fobs.read = read.csv("AAPL forc clean obs.csv", header=TRUE)
fobs = array(unlist(fobs.read$forc.use))

# Reload forecasts
forecast.df = read.csv("AAPL forc_EGIG.TSCD.csv", header=TRUE)

forc1 = forecast.df$x

# ---- Forecast Accuracy Metrics ----
MSFE_1 = mean((fobs - forc1)^2)
MAFE_1 = mean(abs(fobs - forc1))
QLIKE_1 = mean((fobs / forc1) - log(fobs / forc1) - 1, na.rm = TRUE)

Out_1 = rbind(MSFE_1, MAFE_1, QLIKE_1)
colnames(Out_1) = "AAPL EGIG.TSCD"

# Combine into one data frame
Out_combined = data.frame(
  Metric = c("MSFE", "MAFE", "QLIKE"),
  forc1 = Out_1[, 1]
)

# Save to CSV
write.csv(Out_combined, "AAPL Forc Error_EGIG.TSCD.csv", row.names = FALSE)

Out_1



