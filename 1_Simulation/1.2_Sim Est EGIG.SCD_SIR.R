###############################################################################
# This script estimates the EGIG.SCD model parameters using a
# SIR particle filter.
#
# The procedure:
# 1. Reads simulated observation data.
# 2. Approximates the log-likelihood via a particle filter.
# 3. Maximizes the likelihood using L-BFGS-B optimization.
# 4. Repeats estimation for R replications.
# 5. Saves parameter estimates to a CSV file for Monte Carlo analysis.
###############################################################################

rm(list = ls())

# ===================== Read Data =====================
x.read = read.csv("EGIG.SCD Sim_obs_R.csv", header=TRUE)

# ===================== Selected n & N =====================
n = 1000 #sample size (row): 1000, 3000 & 6500
N = 1000 #number of particles (column)
R = ncol(x.read) 

# ===================== Initial Parameter Values =====================
miu.int = 0.6
be.int = 0.97
lam.int = 0.4
del.int = 1.1
w.int = 0.2
sig.int = 0.1
initialTheta = c(miu.int, be.int, lam.int, del.int, w.int, sig.int)
true.pars = c(0.6, 0.97, 0.4, 1.1, 0.2, 0.1)

# ===================== Log-Likelihood via SIR Particle Filter ==================
LL = function(theta, x){  
  miu = theta[1]
  be = theta[2]
  lam = theta[3]
  del = theta[4]
  w = theta[5]
  sig = theta[6]
  
  llike = numeric(n-1)
  
  # Storage objects for particle filter
  un.weight = matrix(data=0, nrow=n-1, ncol=N)
  n.weight = matrix(data=0, nrow=n-1, ncol=N)
  resample.location = matrix(data=0, nrow=n-1, ncol=N)
  resamples = matrix(data=0, nrow=n-1, ncol=N)
  psi = matrix(data=0, nrow=n-1, ncol=N)
  psi.est = matrix(data=0, nrow=n-1, ncol=N)
  
  for(i in 1:(n-1)){
    if(i == 1){
      
      # -------- Draw initial particles from stationary distribution --------
      psi0 = rnorm(n=N, mean=miu, sd=sig/sqrt((1-be^2)))
      
      #Generate state from model:
      psi[i, ] = miu*(1-be) + be*psi0 + rnorm(N, mean=0, sd=sig)
      
    } else {
      
      #Generate state from model:
      psi[i, ] = miu*(1-be) + be*psi.est[i-1, ] + rnorm(N, mean=0, sd=sig)
      
    }
    
    # ================= Likelihood Weight Calculation =================
    #Compute normalized weight:
    b1 = besselK(w, lam/del, expon.scaled = FALSE)
    b2 = besselK(w, (lam+1)/del, expon.scaled = FALSE)
    
    un.weight[i, ] = del/2 * (b2^lam / b1^(lam+1)) *
                    (x.obs[i]^(lam-1) / exp(psi[i, ])^lam) *
                    exp(-w/2 * (((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^del +
                                ((b2 * x.obs[i]) / (b1 * exp(psi[i, ])))^-del))
    
    #Compute normalized weight
    n.weight[i, ] = un.weight[i, ] / sum(un.weight[i, ])
    
    # ================= Resampling Step (SIR) =================
    resample.location[i, ] = sample(x=1:N, replace=TRUE, prob=n.weight[i, ])
    resamples[i, ] = psi[i, resample.location[i, ]]
    
    psi.est[i, ] = resamples[i, ]
    
    # ================= Log-Likelihood Approximation =================
    miu.weigth = 1/N * sum(un.weight[i, ])
    var.weigth = 1/(N-1) * sum((un.weight[i, ] - miu.weigth)^2)
    
    llike[i] = log(miu.weigth) + 1/2 * var.weigth/(N * miu.weigth^2)
  }
  
  LogLike = sum(llike)
  return(-LogLike)
}


# ===================== Optimization Setup =====================
all.out = NULL
outputs = NULL


for(k in 1:R){ 
  
  print(k)
  
  x.obs = tail(array(unlist(x.read[, k])), n)
  
  OP = optimx::optimx(par=initialTheta, fn=LL, gr=NULL, hess=NULL,
                      hessian=FALSE,
                      lower=c(-Inf, -0.99, -Inf, 0, 0, 0),
                      upper=c(Inf, 0.99, Inf, Inf, Inf, Inf),
                      method=c("L-BFGS-B"), 
                      control=list(trace=1, kkt=FALSE))
  
  miu.est = OP$p1
  be.est = OP$p2
  lam.est = OP$p3
  del.est = OP$p4
  w.est = OP$p5
  sig.est = OP$p6
  par.est = c(miu.est, be.est, sig.est, lam.est, del.est, w.est)
  
  all.out = c(par.est)
  outputs = rbind(outputs, all.out)
  
}

# ===================== Save Estimates =====================
write.table(outputs, 
            file="Sim Est_EGIG.SCD_n=1000 (R=1000).csv", 
            # file="Sim Est_EGIG.SCD_n=3000 (R=1000).csv",
            # file="Sim Est_EGIG.SCD_n=6500 (R=1000).csv",
            sep=",", col.names=paste0(c('miu', 'be', 'sig', 'lam', 'del', 'w')), 
            row.names=FALSE)





