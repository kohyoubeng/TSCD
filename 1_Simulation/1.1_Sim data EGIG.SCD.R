###############################################################################
# This script simulates data from the EGIG.SCD model for use in SIR estimation.
#
# 1. Defines model parameters.
# 2. Simulates ψ_t (state), η_t (state error), ε_t (EGIG shock), and y_t (obs).
# 3. Repeats the simulation R times.
# 4. Stores and exports all simulated series into CSV files.
###############################################################################

rm(list = ls())
n = 10000 #Sample size
R = 1000 #Number of replications

# ======================== Model Parameters ========================
miu = 0.6
be = 0.97
lam = 0.4
del = 1.1
w = 0.2
sig = 0.1

# ======================== Storage Matrices ========================
out.obs = matrix(data=NA, nrow=n, ncol=R)
out.eps = matrix(data=NA, nrow=n, ncol=R)
out.eta = matrix(data=NA, nrow=n, ncol=R)
out.psi = matrix(data=NA, nrow=n, ncol=R)

# ===================== Model Simulation Function =====================
fx = function(miu, be, lam, del, w, sig, n){
  
  psi = numeric(n)
  obs = numeric(n)
  
  # -------- Initial state (stationary distribution of AR(1)) --------
  psi0 = rnorm(n=1, mean=miu, sd=sig/sqrt(1-be^2))
  
  # -------- State error η_t ~ N(0, σ²) --------
  eta = rnorm(n=n, mean=0, sd=sig)
  
  # ==================== EGIG Distribution for ε_t ====================
  r = besselK(w, (lam+1)/del, expon.scaled = FALSE)
  s = besselK(w, lam/del, expon.scaled = FALSE)
  v = (r/s)^del 
  
  # -------- Draw EGIG random variables --------
  e.gig = GeneralizedHyperbolic::rgig(n=n, chi=w*v, psi=w/v, lambda=lam/del, 
                                      param=c(chi=w*v, psi=w/v, lambda=lam/del))
  eps = e.gig^(1/del)
  
  # -------- First observation --------
  psi[1] =  miu*(1-be) + be*psi0 + eta[1]
  
  obs[1] = exp(psi[1]) * eps[1]
  
  # -------- Recursive state & observation update --------
  for(i in 2:n){
    
    psi[i] = miu*(1-be) + be*psi[i-1] + eta[i]
    
    obs[i] = exp(psi[i]) * eps[i]
    
  }
  
  my.list = list("obs"=obs, "eps"=eps, "eta"=eta, "psi"=psi)
  return(my.list)
  
}

# ===================== Loop =====================
for(k in 1:R){
  
  Sim.data = fx(miu, be, lam, del, w, sig, n)
  
  out.obs[, k] = Sim.data$obs
  out.eps[, k] = Sim.data$eps
  out.eta[, k] = Sim.data$eta
  out.psi[, k] = Sim.data$psi
  
}

# ===================== Save Results to CSV =====================
write.table(out.obs, file="EGIG.SCD Sim_obs_R.csv", 
            col.names=paste0(rep('obs', R)), sep=",", row.names=FALSE)

write.table(out.eps, file="EGIG.SCD Sim_eps_R.csv", 
            col.names=paste0(rep('eps', R)), sep=",", row.names=FALSE)

write.table(out.eta, file="EGIG.SCD Sim_eta_R.csv", 
            col.names=paste0(rep('eta', R)), sep=",", row.names=FALSE)

write.table(out.psi, file="EGIG.SCD Sim_psi_R.csv", 
            col.names=paste0(rep('psi', R)), sep=",", row.names=FALSE)






