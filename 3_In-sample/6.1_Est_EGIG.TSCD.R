###############################################################################
# REPLICATION INSTRUCTIONS — EGIG.TSCD Model Estimation
#
# This script estimates the EGIG.TSCD duration model using a particle filter
# likelihood and numerical optimization. The estimation procedure, likelihood
# function, and algorithm are identical for AAPL and TSLA. Differences between
# assets arise only from the input data and the choice of initial parameter
# values used to start the optimizer.
#
# --- To reproduce TSLA results ---
# • Replace input file with:
#     "Adj dur TSLA In-sample (t)_29 Sep-10 Oct 2025.csv"
# • Initial parameter values may be adjusted for numerical stability
#
# All other components remain unchanged.
#
###############################################################################
# !!! IMPORTANT !!!
# Sources of script and data **MUST** save in the correct working directory
#
# To obtain replicable results, run: "0.1_Run 500 times.R"
#
# Then use                         : "0.2_Combine all run files.R"
# to combine the outputs from all runs.
###############################################################################

rm(list = setdiff(ls(), c("run", "x.obs", "n", "N")))

# ===================== Load Data =====================
x.read = read.csv("Adj dur AAPL In-sample (t)_29 Sep-10 Oct 2025.csv", 
                  header=TRUE)
x.obs = array(unlist(x.read[, 2]))

n = length(x.obs) # sample size (row)
N = 1000 # number of particles (column)

# ==== Initial Parameter Values (may be adjusted for numerical stability) ======
miu1.int = 0.3
miu2.int = 1
be1.int  = 0.9
be2.int  = 0.95
sig1.int = 0.2
sig2.int = 0.1
lam1.int = 0.05   
lam2.int = 1.1
del1.int = 1.5
del2.int = 0.5
w1.int   = 0.3
w2.int   = 0.1
initialTheta = c(miu1.int, miu2.int, be1.int, be2.int, sig1.int, sig2.int,
                 lam1.int, lam2.int, del1.int, del2.int, w1.int, w2.int)

# ===================== Particle Log-Likelihood Function =====================
LL = function(theta, x){
  miu1 = theta[1]
  miu2 = theta[2]
  be1  = theta[3]
  be2  = theta[4]
  sig1 = theta[5]
  sig2 = theta[6]
  lam1 = theta[7]
  lam2 = theta[8]
  del1 = theta[9]
  del2 = theta[10]
  w1   = theta[11]
  w2   = theta[12]
  
  llike = numeric(n)
  
  un.weight = matrix(0, n, N)
  n.weight  = matrix(0, n, N)
  resample.location = matrix(0, n, N)
  resamples = matrix(0, n, N)
  psi       = matrix(0, n, N)
  psi.est   = matrix(0, n, N)
  
  # Threshold value, r = median
  r = median(x.obs)
  
  for(i in 1:n){
    if(i == 1){
      if(x.obs[i] <= r){
        # -------- Initial latent state from stationary distribution --------
        psi0 = rnorm(N, mean=miu1, sd=sig1/sqrt((1-be1^2)))
        psi[i, ] = miu1*(1-be1) + be1*psi0 + rnorm(N, mean=0, sd=sig1)
        
        # ================= Calculate weight =================
        b1 = besselK(w1, lam1/del1, expon.scaled = FALSE)
        b2 = besselK(w1, (lam1+1)/del1, expon.scaled = FALSE)
        un.weight[i, ] = del1/2 * (b2^lam1 / b1^(lam1+1)) *
          (x.obs[i]^(lam1-1) / exp(psi[i, ])^lam1) *
          exp(-w1/2 * (((b2*x.obs[i])/(b1*exp(psi[i, ])))^del1 +
                         ((b2*x.obs[i])/(b1*exp(psi[i, ])))^-del1))
      } else {
        psi0 = rnorm(N, mean=miu2, sd=sig2/sqrt((1-be2^2)))
        psi[i, ] = miu2*(1-be2) + be2*psi0 + rnorm(N, mean=0, sd=sig2)
        
        b1 = besselK(w2, lam2/del2, expon.scaled = FALSE)
        b2 = besselK(w2, (lam2+1)/del2, expon.scaled = FALSE)
        un.weight[i, ] = del2/2 * (b2^lam2 / b1^(lam2+1)) *
          (x.obs[i]^(lam2-1) / exp(psi[i, ])^lam2) *
          exp(-w2/2 * (((b2*x.obs[i])/(b1*exp(psi[i, ])))^del2 +
                         ((b2*x.obs[i])/(b1*exp(psi[i, ])))^-del2))
      }
    } else {
      if(x.obs[i-1] <= r){
        psi[i, ] = miu1*(1-be1) + be1*psi.est[i-1, ] + rnorm(N, mean=0, sd=sig1)
        
        b1 = besselK(w1, lam1/del1, expon.scaled = FALSE)
        b2 = besselK(w1, (lam1+1)/del1, expon.scaled = FALSE)
        un.weight[i, ] = del1/2 * (b2^lam1 / b1^(lam1+1)) *
          (x.obs[i]^(lam1-1) / exp(psi[i, ])^lam1) *
          exp(-w1/2 * (((b2*x.obs[i])/(b1*exp(psi[i, ])))^del1 +
                         ((b2*x.obs[i])/(b1*exp(psi[i, ])))^-del1))
      } else {
        psi[i, ] = miu2*(1-be2) + be2*psi.est[i-1, ] + rnorm(N, mean=0, sd=sig2)
        
        b1 = besselK(w2, lam2/del2, expon.scaled = FALSE)
        b2 = besselK(w2, (lam2+1)/del2, expon.scaled = FALSE)
        un.weight[i, ] = del2/2 * (b2^lam2 / b1^(lam2+1)) *
          (x.obs[i]^(lam2-1) / exp(psi[i, ])^lam2) *
          exp(-w2/2 * (((b2*x.obs[i])/(b1*exp(psi[i, ])))^del2 +
                         ((b2*x.obs[i])/(b1*exp(psi[i, ])))^-del2))
      }
    }
    
    # Normalised weight
    n.weight[i, ] = un.weight[i, ] / sum(un.weight[i, ])
    
    # ================= Resampling (SIR) =================
    resample.location[i, ] = sample(1:N, size = N, replace = TRUE, 
                                    prob = n.weight[i, ])
    psi.est[i, ] = psi[i, resample.location[i, ]]
    
    # ================= Log-Likelihood Approximation =================
    miu.weigth = mean(un.weight[i, ])
    var.weigth = var(un.weight[i, ])
    
    llike[i] = log(miu.weigth) + 0.5 * var.weigth/(N * miu.weigth^2)
  }
  
  LogLike = sum(llike)
  return(-LogLike)
}

# ===================== Optimization =====================
OP = optimx::optimx(par=initialTheta,
                    fn=LL,
                    gr=NULL,
                    hess=NULL,
                    hessian=FALSE,
                    lower=c(-Inf, -Inf, -0.99, -0.99, 0, 0, -Inf, -Inf, 0, 0, 
                            0, 0),
                    upper=c( Inf,  Inf,  0.99,  0.99, Inf, Inf, Inf, Inf, Inf, 
                             Inf, Inf, Inf),
                    method=c("L-BFGS-B"),
                    control=list(trace=1, kkt=FALSE))

theta_hat_raw = as.numeric(OP[1, 1:length(initialTheta)])

# ===================== Model Selection Criteria =====================
param_names = c("miu1", "miu2", "be1", "be2", "sig1", "sig2",
                "lam1", "lam2", "del1", "del2", "w1", "w2")

LogL = -OP$value
aic  = 2*length(initialTheta) - 2*LogL
bic  = length(initialTheta) * log(n) - 2*LogL

# ===================== Output Results =====================
results = data.frame(
  Parameter = "Estimate",
  miu1 = theta_hat_raw[1],
  miu2 = theta_hat_raw[2],
  be1  = theta_hat_raw[3],
  be2  = theta_hat_raw[4],
  sig1 = theta_hat_raw[5],
  sig2 = theta_hat_raw[6],
  lam1 = theta_hat_raw[7],
  lam2 = theta_hat_raw[8],
  del1 = theta_hat_raw[9],
  del2 = theta_hat_raw[10],
  w1   = theta_hat_raw[11],
  w2   = theta_hat_raw[12],
  LogL = LogL,
  AIC  = aic,
  BIC  = bic,
  stringsAsFactors = FALSE
)

print(results)

