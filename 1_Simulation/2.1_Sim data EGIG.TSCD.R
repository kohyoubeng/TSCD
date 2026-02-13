###############################################################################
# This script simulates data from the EGIG.TSCD model for use in SIR estimation.
#
# 1. Defines model parameters.
# 2. Simulates ψ_t (state), η_t (state error), ε_t (EGIG shock), and y_t (obs).
# 3. Applies regime switching based on a threshold rule.
# 4. Repeats the simulation R times.
# 5. Stores and exports all simulated series into CSV files.
###############################################################################

rm(list = ls())
n = 10000 #Sample size
R = 1000 #Number of replications

# ======================== Model Parameters ========================
miu1 = 0.3
miu2 = 1.2
be1 = 0.9
be2 = 0.95
lam1 = 0.3
lam2 = 0.5
del1 = 1.05
del2 = 1.3
w1 = 0.3
w2 = 0.1
sig1 = 0.15
sig2 = 0.06

# ======================== Storage Matrices ========================
out.obs = matrix(data=NA, nrow=n, ncol=R)
out.eps = matrix(data=NA, nrow=n, ncol=R)
out.eta = matrix(data=NA, nrow=n, ncol=R)
out.psi = matrix(data=NA, nrow=n, ncol=R)

# ===================== Model Simulation Function =====================
fx = function(miu1, miu2, be1, be2, lam1, lam2, del1, del2, w1, w2, sig1, sig2, 
              n){
  
  eta = numeric(n)
  eps = numeric(n)
  psi = numeric(n)
  obs = numeric(n)
  
  #Let initial obs, obs0=1 and threshold value, r=1.12
  obs0 = 1
  r = 1.12
  
  for(i in 1:n){
    if(i == 1){
      
      # -------- Initial state (stationary distribution of AR(1)) --------
      psi0 = rnorm(n=1, mean=miu1, sd=sig1/sqrt(1-be1^2))
      
      # -------- State error η_t ~ N(0, σ²) --------
      eta[1] = rnorm(n=1, mean=0, sd=sig1)
      
      # ==================== EGIG Distribution for ε_t ====================
      p = besselK(w1, lam1/del1, expon.scaled = FALSE)
      q = besselK(w1, (lam1+1)/del1, expon.scaled = FALSE)
      v = (p/q)^del1
      
      # -------- Draw EGIG random variables --------
      e.gig = GeneralizedHyperbolic::rgig(n=1, chi=w1*v, psi=w1/v, 
                                  lambda=lam1/del1, 
                                  param=c(chi=w1*v, psi=w1/v, lambda=lam1/del1)) 
      eps[1] = e.gig^(1/del1)
      
      
      # -------- First observation --------
      psi[1] =  miu1*(1-be1) + be1*psi0 + eta[1]
      
      obs[1] = exp(psi[1]) * eps[1]
      
    } else {
      if(obs[i-1] <= r){
        
        #State error, eta~N(0, sig^2)
        eta[i] = rnorm(n=1, mean=0, sd=sig1)
        
        #EGIG distribution
        p = besselK(w1, lam1/del1, expon.scaled = FALSE)
        q = besselK(w1, (lam1+1)/del1, expon.scaled = FALSE)
        v = (p/q)^del1
        
        e.gig = GeneralizedHyperbolic::rgig(n=1, chi=w1*v, psi=w1/v, 
                                  lambda=lam1/del1, 
                                  param=c(chi=w1*v, psi=w1/v, lambda=lam1/del1)) 
        
        eps[i] = e.gig^(1/del1)
        
        psi[i] =  miu1*(1-be1) + be1*psi[i-1] + eta[i]
        
        obs[i] = exp(psi[i]) * eps[i]
        
      } else {
        
        #State error, eta~N(0, sig^2)
        eta[i] = rnorm(n=1, mean=0, sd=sig2)
        
        #EGIG distribution
        p = besselK(w2, lam2/del2, expon.scaled = FALSE)
        q = besselK(w2, (lam2+1)/del2, expon.scaled = FALSE)
        v = (p/q)^del2
        
        e.gig = GeneralizedHyperbolic::rgig(n=1, chi=w2*v, psi=w2/v, 
                                  lambda=lam2/del2, 
                                  param=c(chi=w2*v, psi=w2/v, lambda=lam2/del2)) 
        eps[i] = e.gig^(1/del2)
        
        psi[i] =  miu2*(1-be2) + be2*psi[i-1] + eta[i]
  
        obs[i] = exp(psi[i]) * eps[i]
        
      }
    }
  }
  
  my.list = list("obs"=obs, "eps"=eps, "eta"=eta, "psi"=psi)
  return(my.list)
  
}

# ===================== Loop =====================
for(k in 1:R){
  
  print(k)
  
  Sim.data = fx(miu1, miu2, be1, be2, lam1, lam2, del1, del2, w1, w2, 
                sig1, sig2, n)
  
  out.obs[, k] = Sim.data$obs
  out.eps[, k] = Sim.data$eps
  out.eta[, k] = Sim.data$eta
  out.psi[, k] = Sim.data$psi
  
}

# ===================== Save Results to CSV =====================
write.table(out.obs, file="EGIG.TSCD Sim_obs_R.csv", 
            col.names=paste0(rep('obs', R)), sep=",", row.names=FALSE)

write.table(out.eps, file="EGIG.TSCD Sim_eps_R.csv", 
            col.names=paste0(rep('eps', R)), sep=",", row.names=FALSE)

write.table(out.eta, file="EGIG.TSCD Sim_eta_R.csv", 
            col.names=paste0(rep('eta', R)), sep=",", row.names=FALSE)

write.table(out.psi, file="EGIG.TSCD Sim_psi_R.csv", 
            col.names=paste0(rep('psi', R)), sep=",", row.names=FALSE)





