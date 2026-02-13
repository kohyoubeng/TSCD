###############################################################################
# This script estimates the EGIG.TSCD model parameters using a
# SIR particle filter.
#
# The procedure:
# 1. Reads simulated observation data.
# 2. Approximates the log-likelihood via a particle filter with
#    regime switching based on a threshold rule.
# 3. Maximizes the likelihood using L-BFGS-B optimization.
# 4. Repeats estimation for R replications.
# 5. Saves parameter estimates to a CSV file for Monte Carlo analysis.
###############################################################################

rm(list = ls())

# ===================== Read Data =====================
x.read = read.csv("EGIG.TSCD Sim_obs_R.csv", header=TRUE)

# ===================== Selected n & N =====================
n = 3000 #sample size (row): 1000, 3000 & 6500
N = 1000 #number of particles (column)
R = ncol(x.read) 

# ===================== Initial Parameter Values =====================
miu1.int = 0.3	
miu2.int = 1.1
be1.int = 0.9
be2.int = 0.95	
lam1.int = 0.3
lam2.int = 0.5	
del1.int = 1.05
del2.int = 1.3
w1.int = 0.3
w2.int = 0.1	
sig1.int = 0.15
sig2.int = 0.06
initialTheta = c(miu1.int, miu2.int, be1.int, be2.int, lam1.int, lam2.int, 
                 del1.int, del2.int, w1.int, w2.int, sig1.int, sig2.int)

# ===================== Log-Likelihood via SIR Particle Filter ==================
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
  
  llike = numeric(n-1)
  
  # Storage objects for particle filter
  un.weight = matrix(data=0, nrow=n-1, ncol=N)
  n.weight = matrix(data=0, nrow=n-1, ncol=N)
  resample.location = matrix(data=0, nrow=n-1, ncol=N)
  resamples = matrix(data=0, nrow=n-1, ncol=N)
  psi = matrix(data=0, nrow=n-1, ncol=N)
  psi.est = matrix(data=0, nrow=n-1, ncol=N)
  
  #threshold value, r=median
  r = median(x.obs)
  
  for(i in 1:(n-1)){
    if(i == 1){
      if(x.obs[i] <= r){
        
        # -------- Draw initial particles from stationary distribution --------
        #Give samples:
        psi0 = rnorm(n=N, mean=miu1, sd=sig1/sqrt((1-be1^2)))
        
        #Generate state from model:
        psi[i, ] = miu1*(1-be1) + be1*psi0 + rnorm(N, mean=0, sd=sig1)
        
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
        
        # ================= Likelihood Weight Calculation =================
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
        
        # ================= Likelihood Weight Calculation =================
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
        
        # ================= Likelihood Weight Calculation =================
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
          lower=c(-Inf, -Inf, -0.99, -0.99, -Inf, -Inf, 0, 0, 0, 0, 0, 0),
          upper=c(Inf, Inf, 0.99, 0.99, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
          method=c("L-BFGS-B"), control=list(trace=1, kkt=FALSE))
  
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
  
  all.out = c(par.est)
  outputs = rbind(outputs, all.out)
  
}

# ===================== Save Estimates =====================
write.table(outputs, 
            # file="Sim Est_EGIG.TSCD_n=1000 (R=1000).csv", 
            file="Sim Est_EGIG.TSCD_n=3000 (R=1000).csv",
            # file="Sim Est_EGIG.TSCD_n=6500 (R=1000).csv",
            sep=",",
            col.names=paste0(c('miu1', 'miu2', 'be1', 'be2', 'sig1', 'sig2', 
                               'lam1', 'lam2', 'del1', 'del2', 'w1', 'w2')), 
            row.names=FALSE)
