###############################################################################
# PSI ESTIMATION — Compute latent states for Burr.SCD model
#
# This script applies a particle filter to compute the latent states (ψ_t) 
# and their exponentiated values (exp(ψ_t)) for a chosen stock (AAPL or TSLA):
#   • Reads in in-sample observed data
#   • Uses parameter values obtained from the estimation step 
#     (e.g., miu, be, sig, a, q)
#   • Applies particle filtering with resampling to estimate ψ_t
#   • Computes exp(ψ_t) for each time point
#   • Saves the results as a CSV file combining ψ_t and exp(ψ_t)
#
# --- Notes ---
# • MUST get the output CSV files of estimated parameter values for AAPL or TSLA
# • Change the input CSV files and parameter values for AAPL or TSLA
# • Output files are named like "AAPL psi & exp(psi)_Burr.SCD.csv" 
#   or "TSLA psi & exp(psi)_Burr.SCD.csv"
###############################################################################

rm(list = ls())

# ===================== Read Data =====================
x.read = read.csv("Adj dur AAPL In-sample (t)_29 Sep-10 Oct 2025.csv", 
                  header=TRUE)
par.read = read.csv("AAPL Est Burr.SCD_500.csv", header=TRUE)

x.obs = array(unlist(x.read[, 2]))

n = nrow(x.read) #sample size (row)
N = 1000         #number of particles (column)

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
my_func = function(x.obs, miu, be, sig, a, q){
  
  for(i in 1:n){
    if(i == 1){
      
      # -------- Initial latent state from stationary AR(1) distribution --------
      psi0 = rnorm(n=N, mean=miu, sd=sig/sqrt((1-be^2)))
      
      #Generate state from model: generate 1st state
      psi[i, ] = miu*(1-be) + be*psi0 + rnorm(N, mean=0, sd=sig)
      
    } else {
      
      #Generate state from model: generate 2nd state & so on
      psi[i, ] = miu*(1-be) + be*psi.est[i-1, ] + rnorm(N, mean=0, sd=sig)
      
    }
    
    # ================= Observation Density =================
    #Compute normalized weight:
    B = beta(1+(1/a), q-(1/a))
    # B = beta(1+(1/a), b)
    
    un.weight[i, ] = a * q^(a+1) * B^a *
      (1 / (1 + (q * B * x.obs[i] / exp(psi[i, ]))^a)^(q+1)) *
      x.obs[i]^(a-1) * (1/exp(psi[i, ]))^a
    
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
miu = 0.6724
be = 0.9095
sig = 0.2489
a = 1.3567
q = 2.1860

psi = my_func(x.obs, miu, be, sig, a, q)
phi = exp(psi)

out = cbind(psi, phi)

write.csv(out, "AAPL psi & exp(psi)_Burr.SCD.csv", row.names=F)



