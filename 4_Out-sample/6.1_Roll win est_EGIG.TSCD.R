###############################################################################
# Rolling Window Parameter Estimation — EGIG TSCD Model (Particle Filter, SIR)
#
# This script performs rolling-window maximum likelihood estimation of the
# EGIG–TSCD model using a SIR particle filter to approximate the likelihood.
#
# The procedure:
# • Uses a particle filter (N = 1000) to evaluate the likelihood
# • Estimates parameters (μ, β, σ, a) in each rolling window
# • Maximizes the particle-approximated log-likelihood via L-BFGS-B
# • Saves parameter estimates and log-likelihood for each window separately
#
# The SAME script applies to both AAPL and TSLA datasets.
#
# --- To use for TSLA ---
# • Replace input file:
#     "AAPL x.use roll.csv"  →  "TSLA x.use roll.csv"
# • Change initial parameter vector `init1` to the TSLA starting values
# • Replace output filenames beginning with "AAPL" to "TSLA"
#
# Notes:
# • MUST use rolling-window standardized durations (x.use) as input
# • Each column represents one estimation window
# • Last column is excluded because it has no next-day forecast target
# • Particle filter settings, model structure, and estimation method remain identical
###############################################################################

rm(list = ls())

# ===================== Load Data =====================
x.read = read.csv("AAPL x.use roll.csv", header=TRUE)

n = nrow(x.read)        #sample size (row)
N = 1000                #number of particles (column)

n.f = ncol(x.read) -1 #because the last column not use to forc next day 

# ===================== Particle Log-Likelihood Function =====================
#Loglike function
LL = function(theta, x){
  miu1 = theta[1]
  miu2 = theta[2]
  be1 = theta[3]
  be2 = theta[4]
  lam1 = theta[5]
  lam2 = theta[6]
  del1 = theta[7]
  del2 = theta[8]
  w1 = theta[9]
  w2 = theta[10]
  sig1 = theta[11]
  sig2 = theta[12]
  
  #for all
  llike = numeric(n)
  
  un.weight = matrix(data=0, nrow=n, ncol=N)
  n.weight = matrix(data=0, nrow=n, ncol=N)
  resample.location = matrix(data=0, nrow=n, ncol=N)
  resamples = matrix(data=0, nrow=n, ncol=N)
  psi = matrix(data=0, nrow=n, ncol=N)
  psi.est = matrix(data=0, nrow=n, ncol=N)
  
  #threshold value, r=median
  r = median(x.obs)
  
  for(i in 1:n){
    if(i == 1){
      if(x.obs[i] <= r){
        
        # -------- Initial latent state from stationary distribution --------
        psi0 = rnorm(n=N, mean=miu1, sd=sig1/sqrt((1-be1^2)))
        
        #Generate state from model:
        psi[i, ] = miu1*(1-be1) + be1*psi0 + rnorm(N, mean=0, sd=sig1)
        
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
        
        #Generate state from model:
        psi[i, ] = miu1*(1-be1) + be1*psi.est[i-1, ] + rnorm(N, mean=0, sd=sig1)
        
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
    
    # ================= Resampling (SIR) =================
    resample.location[i, ] = sample(1:N, size = N, replace = TRUE, 
                                    prob = n.weight[i, ])
    psi.est[i, ] = psi[i, resample.location[i, ]]
    
    # ================= Log-Likelihood Approximation =================
    miu.weigth = 1/N * sum(un.weight[i, ])
    var.weigth = 1/(N-1) * sum((un.weight[i, ] - miu.weigth)^2)
    
    llike[i] = log(miu.weigth) + 1/2 * var.weigth/(N * miu.weigth^2)
  }
  
  LogLike = sum(llike)
  return(-LogLike)
}

# =================== Storage for rolling estimation results ===================
out.par.est = NULL
out.LL = NULL
out.all = NULL

# ================= Rolling window estimation =================
for(k in 1:n.f){ 
  
  print(paste("Estimating window", k, "of", n.f))
  
  #Read data
  x.obs = x.read[, k]
  
  # Initial for optimizer
  #miu1,	miu2,	be1,	be2, lam1,	lam2,	del1,	del2,	w1,	w2,	sig1,	sig2	
  init1 = c(0.3, 1.2, 0.93, 0.97, 0.3, 0.5, 1.05, 1.3, 0.3, 0.1, 0.15, 0.06)
  
  # Parameter estimation via maximum likelihood
  OP = optimx::optimx(par=init1, fn=LL, gr=NULL, hess=NULL, hessian=FALSE,
                      lower=c(-Inf, -Inf, -0.99, -0.99, -Inf, -Inf, 0, 0, 0, 0, 0, 0),
                      upper=c(Inf, Inf, 0.99, 0.99, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
                      method=c("L-BFGS-B"), control=list(trace=1, kkt=FALSE))
  
  # Extract parameter estimates
  miu1.est = OP$p1
  miu2.est = OP$p2
  be1.est = OP$p3
  be2.est = OP$p4
  lam1.est = OP$p5
  lam2.est = OP$p6
  del1.est = OP$p7
  del2.est = OP$p8
  w1.est = OP$p9
  w2.est = OP$p10
  sig1.est = OP$p11
  sig2.est = OP$p12
  par.est = c(miu1.est, miu2.est, be1.est, be2.est, sig1.est, sig2.est, 
              lam1.est, lam2.est, del1.est, del2.est, w1.est, w2.est)
  
  # Extract log-likelihood
  LogL = -OP$value
  
  # Store cumulative results
  out.par.est = rbind(out.par.est, par.est)
  out.LL = rbind(out.LL, LogL)
  out.all = cbind(out.par.est, out.LL)
  
  #------------------------------------------
  # ✅ Save each result immediately
  #------------------------------------------
  file.name = paste0("AAPL roll pars_EGIG.TSCD_win", k, ".csv")
  single.out = data.frame(
    window = k,
    miu1 = miu1.est,
    miu2 = miu2.est,
    be1 = be1.est,
    be2 = be2.est,
    sig1 = sig1.est,
    sig2 = sig2.est,
    lam1 = lam1.est,
    lam2 = lam2.est,
    del1 = del1.est,
    del2 = del2.est,
    w1 = w1.est,
    w2 = w2.est,
    LogL = LogL
  )
  write.csv(single.out, file.name, row.names = FALSE)
  #------------------------------------------
  
}

