###############################################################################
# Rolling Window Parameter Estimation — Wei SCD Model (Particle Filter, SIR)
#
# This script performs rolling-window maximum likelihood estimation of the
# Wei–SCD model using a SIR particle filter to approximate the likelihood.
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
LL = function(theta){  
  miu = theta[1]
  be = theta[2]
  a = theta[3]
  sig = theta[4]
  
  #for all
  llike = numeric(n)
  
  un.weight = matrix(data=0, nrow=n, ncol=N)
  n.weight = matrix(data=0, nrow=n, ncol=N)
  resample.location = matrix(data=0, nrow=n, ncol=N)
  resamples = matrix(data=0, nrow=n, ncol=N)
  psi = matrix(data=0, nrow=n, ncol=N)
  psi.est = matrix(data=0, nrow=n, ncol=N)
  
  for(i in 1:n) {
    
    if(i == 1){
      
      # -------- Initial latent state from stationary distribution --------
      psi0 = rnorm(n=N, mean=miu, sd = sig/sqrt(1-be^2))
      
      #Generate state from model:
      psi[i, ] = miu*(1-be) + be*psi0 + rnorm(N, mean=0, sd=sig)
      
    } else {
      
      #Generate state from model:
      psi[i, ] = miu*(1-be) + be*psi.est[i-1, ] + rnorm(N, mean=0, sd=sig)
      
    }
    
    # ================= Calculate weight =================
    #Compute unnormalized weight:
    g = gamma(1 + 1/a)
    
    un.weight[i, ] = a * g^a * x.obs[i]^(a - 1) / (exp(psi[i, ]))^a * 
      exp(-((x.obs[i] * g) / exp(psi[i, ]))^a)
    
    #Compute normalized weight:
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
  
  # Read data
  x.obs = x.read[, k]
  
  # Initial for optimizer
  init1 = c(0.6, 0.9, 0.5, 0.3)
  
  # Parameter estimation via maximum likelihood
  OP = optimx::optimx(par=init1, fn=LL, gr=NULL, hess=NULL, hessian=FALSE,
                      lower=c(-Inf, -0.99, 0, 0),
                      upper=c(Inf, 0.99, Inf, Inf),
                      method=c("L-BFGS-B"),
                      control=list(trace=0, kkt=FALSE))
  
  # Extract parameter estimates
  miu.est = OP$p1
  be.est  = OP$p2
  a.est   = OP$p3
  sig.est = OP$p4
  par.est = c(miu.est, be.est, sig.est, a.est)
  
  # Extract log-likelihood
  LogL = -OP$value
  
  # Store cumulative results
  out.par.est = rbind(out.par.est, par.est)
  out.LL = rbind(out.LL, LogL)
  out.all = cbind(out.par.est, out.LL)
  
  #------------------------------------------
  # ✅ Save each result immediately
  #------------------------------------------
  file.name = paste0("AAPL roll pars_Wei.SCD_win", k, ".csv")
  single.out = data.frame(
    window = k,
    miu = miu.est,
    be = be.est,
    sig = sig.est,
    a = a.est,
    LogL = LogL
  )
  write.csv(single.out, file.name, row.names = FALSE)
  #------------------------------------------
}

