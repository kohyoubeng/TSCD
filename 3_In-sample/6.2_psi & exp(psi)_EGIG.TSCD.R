###############################################################################
# PSI ESTIMATION — Compute latent states for EGIG.TSCD model
#
# This script applies a particle filter to compute the latent states (ψ_t) 
# and their exponentiated values (exp(ψ_t)) for a chosen stock (AAPL or TSLA):
#   • Reads in in-sample observed data
#   • Uses parameter values obtained from the estimation step 
#     (e.g., r, miu1, miu2, be1, be2, sig1, sig2, lam1, lam2, del1, del2, w1, w2)
#   • Applies particle filtering with resampling to estimate ψ_t
#   • Computes exp(ψ_t) for each time point
#   • Saves the results as a CSV file combining ψ_t and exp(ψ_t)
#
# --- Notes ---
# • MUST get the output CSV files of estimated parameter values for AAPL or TSLA
# • Change the input CSV files and parameter values for AAPL or TSLA
# • Output files are named like "AAPL psi & exp(psi)_EGIG.TSCD.csv" 
#   or "TSLA psi & exp(psi)_EGIG.TSCD.csv"
###############################################################################

rm(list = ls())

# ===================== Read Data =====================
x.read = read.csv("Adj dur AAPL In-sample (t)_29 Sep-10 Oct 2025.csv", 
                  header=TRUE)
par.read = read.csv("AAPL Est EGIG.TSCD_500.csv", header=TRUE)

x.obs = array(unlist(x.read[, 2]))

n = nrow(x.read) #sample size (row)
N = 1000         #number of particles (column)
r = median(x.obs)

# ===================== Particle Filter Storage =====================
un.weight = matrix(data=0, nrow=n, ncol=N)
n.weight = matrix(data=0, nrow=n, ncol=N)
resample.location = matrix(data=0, nrow=n, ncol=N)
resamples = matrix(data=0, nrow=n, ncol=N)
psi = matrix(data=0, nrow=n, ncol=N)
psi.est = matrix(data=0, nrow=n, ncol=N)

epsi = numeric(n)

out.epsi = NULL

# ===================== Particle Filtering Function =====================
my_func = function(x.obs, r, miu1, miu2, be1, be2, sig1, sig2, lam1, lam2, 
                   del1, del2, w1, w2){
  
  for(i in 1:n){
    if(i == 1){
      if(x.obs[i] <= r){
        
        # -------- Initial latent state from stationary AR(1) distribution --------
        psi0 = rnorm(n=N, mean=miu1, sd=sig1/sqrt((1-be1^2)))
        
        #Generate state from model: generate 1st state
        psi[i, ] = miu1*(1-be1) + be1*psi0 + rnorm(N, mean=0, sd=sig1)
        
        # ================= Observation Density =================
        #Compute normalized weight:
        b1 = besselK(w1, lam1/del1, expon.scaled = FALSE)
        b2 = besselK(w1, (lam1+1)/del1, expon.scaled = FALSE)
        
        un.weight[i, ] = del1/2 * (b2^lam1 / b1^(lam1+1)) *
          (x.obs[i]^(lam1-1) / exp(psi[i, ])^lam1) *
          exp(-w1/2 * (((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^del1 +
                         ((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^-del1))
        
      } else {
        
        #Give samples:
        psi0 = rnorm(n=N, mean=miu2, sd=sig2/sqrt((1-be2^2)))
        
        #Generate state from model:
        psi[i, ] = miu2*(1-be2) + be2*psi0 + rnorm(N, mean=0, sd=sig2)
        
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
        
        #Generate state from model: generate 2nd state & so on
        psi[i, ] = miu1*(1-be1) + be1*psi.est[i-1, ] + rnorm(N, mean=0, sd=sig1)
        
        # ================= Observation Density =================
        #Compute normalized weight:
        b1 = besselK(w1, lam1/del1, expon.scaled = FALSE)
        b2 = besselK(w1, (lam1+1)/del1, expon.scaled = FALSE)
        
        un.weight[i, ] = del1/2 * (b2^lam1 / b1^(lam1+1)) *
          (x.obs[i]^(lam1-1) / exp(psi[i, ])^lam1) *
          exp(-w1/2 * (((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^del1 +
                         ((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^-del1))
        
      } else {
        
        #Generate state from model:
        psi[i, ] = miu2*(1-be2) + be2*psi.est[i-1, ] + rnorm(N, mean=0, sd=sig2)
        
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
    n.weight[i, ] = un.weight[i, ] / sum(un.weight[i, ])
    
    # ================= Resampling Step =================
    resample.location[i, ] = sample(x=1:N, replace=TRUE, prob=n.weight[i, ])
    resamples[i, ] = psi[i, resample.location[i, ]]
    
    psi.est[i, ] = resamples[i, ]
    
    # ================= State Estimate =================
    epsi[i] = psi[i, ]  %*% n.weight[i, ]
    
  }
  
  return(epsi)
  
}

# ===================== Parameter Values (from estimation) =====================
miu1 = 0.3275
miu2 = 0.9883
be1 = 0.8986
be2 = 0.9571
sig1 = 0.2297
sig2 = 0.1477
lam1 = 0.0542
lam2 = 1.0754
del1 = 1.8206
del2 = 0.4572
w1 = 0.3249
w2 = 0.0934

psi = my_func(x.obs, r, miu1, miu2, be1, be2, sig1, sig2, lam1, lam2, 
              del1, del2, w1, w2)
phi = exp(psi)

out = cbind(psi, phi)

write.csv(out, "AAPL psi & exp(psi)_EGIG.TSCD.csv", row.names=F)



