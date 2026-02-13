###############################################################################
# REPLICATION INSTRUCTIONS — EGIG.SCD Model Estimation
#
# This script estimates the EGIG.SCD duration model using a particle filter
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
miu.int = 0.6
be.int  = 0.97
lam.int = 0.01   
del.int = 1.1
w.int   = 0.2
sig.int = 0.1
initialTheta = c(miu.int, be.int, sig.int, lam.int, del.int, w.int)

# ===================== Particle Log-Likelihood Function =====================
LL = function(theta, x){  
  miu = theta[1]
  be  = theta[2]
  sig = theta[3]
  lam = theta[4]
  del = theta[5]
  w   = theta[6]
  
  llike = numeric(n)
  
  un.weight = matrix(0, nrow=n, ncol=N)
  n.weight  = matrix(0, nrow=n, ncol=N)
  resample.location = matrix(0, nrow=n, ncol=N)
  resamples = matrix(0, nrow=n, ncol=N)
  psi = matrix(0, nrow=n, ncol=N)
  psi.est = matrix(0, nrow=n, ncol=N)
  
  for(i in 1:n){
    if(i == 1){
      # -------- Initial latent state from stationary distribution --------
      psi0 = rnorm(n=N, mean=miu, sd=sig/sqrt((1-be^2)))
      psi[i, ] = miu*(1-be) + be*psi0 + rnorm(N, mean=0, sd=sig)
    } else {
      psi[i, ] = miu*(1-be) + be*psi.est[i-1, ] + rnorm(N, mean=0, sd=sig)
    }
    
    # ================= Calculate weight =================
    b1 = besselK(w, lam/del, expon.scaled = FALSE)
    b2 = besselK(w, (lam+1)/del, expon.scaled = FALSE)
    
    un.weight[i, ] = del/2 * (b2^lam / b1^(lam+1)) *
      (x[i]^(lam-1) / exp(psi[i, ])^lam) *
      exp(-w/2 * (((b2 * x[i]) / (b1 * exp(psi[i, ])))^del +
                    ((b2 * x[i]) / (b1 * exp(psi[i, ])))^-del))
    
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
                    fn=function(th) LL(th, x.obs),
                    gr=NULL,
                    hess=NULL,
                    hessian=FALSE,
                    lower=c(-Inf, -0.99, 0, -Inf, 0, 0),
                    upper=c(Inf, 0.99, Inf, Inf, Inf, Inf),
                    method=c("L-BFGS-B"), 
                    control=list(trace=1, kkt=FALSE))

theta_hat_raw = as.numeric(OP[1, 1:length(initialTheta)])

# ===================== Model Selection Criteria =====================
param_names = c("miu", "be", "sig", "lam", "del", "w")
LogL        = -OP$value
aic         = 2*length(initialTheta) - 2*LogL
bic         = length(initialTheta) * log(n) - 2*LogL

# ===================== Output Results =====================
results = data.frame(Parameter = param_names,
                     Estimate  = theta_hat_raw,
                     LogL      = LogL,
                     AIC       = aic,
                     BIC       = bic,
                     stringsAsFactors = FALSE)

results = data.frame(
  Parameter = "Estimate",
  miu  = theta_hat_raw[1],
  be   = theta_hat_raw[2],
  sig  = theta_hat_raw[3],
  lam  = theta_hat_raw[4],
  del  = theta_hat_raw[5],
  w    = theta_hat_raw[6],
  LogL = LogL,
  AIC  = aic,
  BIC  = bic,
  stringsAsFactors = FALSE
)

print(results)


