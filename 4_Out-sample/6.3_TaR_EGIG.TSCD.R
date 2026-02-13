###############################################################################
# TaR Forecasting and Backtesting — EGIG TSCD Model
#
# This script produces TaR duration forecasts and evaluates
# their reliability using violation analysis and statistical backtests.
#
# The procedure:
# • Computes TaR forecasts at multiple upper risk levels
#   (50%, 75%, 90%, 95%, 97.5%, 99%)
# • Uses forecasted conditional scale (exp(ψ)) and estimated shape parameter a
# • Calculates Violation Rate (VR) and Ratio of Violation Rate (RVR)
# • Performs three standard risk backtests:
#     – Kupiec Unconditional Coverage Test
#     – Dynamic Quantile (DQ) Test
#     – Christoffersen Conditional Coverage Test
# • Saves TaR forecasts, violation statistics, and test results to CSV files
#
# The SAME script applies to both AAPL and TSLA datasets.
#
# --- To use for TSLA ---
# • Replace input file:
#     "AAPL forc.use, SCDs, TSCD & r.csv"  →  "TSLA forc.use, SCDs, TSCD & r.csv"
# • Replace parameter file:
#     "AAPL Roll pars_EGIG.TSCD.csv" → "TSLA Roll pars_EGIG.TSCD.csv"
# • Replace output filenames beginning with "AAPL" to "TSLA"
#
# Notes:
# • MUST include the following columns in the input file:
#     forc.use, Wei, GG, Burr, GB2, EGIG, TSCD, r
# • Wei, GG, Burr, GB2, EGIG & TSCD columns = forecasted conditional scale exp(ψ)
###############################################################################

rm(list = ls())

#TaR risk level  :  50%, 75%,  90%, 95%,  97.5% & 99%
#upper risk level:  0.5, 0.25, 0.1, 0.05, 0.025 & 0.01
u = c(0.5, 0.75, 0.9, 0.95, 0.975, 0.99) #as percentile in calculation

out.fTaR = NULL
out.forc = NULL

#read data: forc.use, exp(fpsi) models & r
#col headings: forc.use, Wei, GG, Burr, GB2, EGIG, TSCD, r
data = read.csv("AAPL forc.use, SCDs, TSCD & r.csv", header=TRUE)

obs.dur = data$forc.use
fphi = data$TSCD  #fphi = exp(fpsi)
r = data$r

#read forecasted parameters
fpar.read = read.csv("AAPL Roll pars_EGIG.TSCD.csv", header=TRUE)

n.f = nrow(fpar.read)
fTaR = numeric(n.f)

###############################################################################
#functions to find z in TaR ##################################################
findprob = function(f, interval, target) {
  optimize(function(y1) {
    abs(integrate(f, 0, y1)$value-target)
  }, interval)$minimum
}

myfunc1 = function(y2){
  y2^(lam1/del1 - 1) * exp(-y2 - w1^2/4 * y2^(-1))
}

myfunc2 = function(y2){
  y2^(lam2/del2 - 1) * exp(-y2 - w2^2/4 * y2^(-1))
}
###############################################################################

# ================= TaR Forecast Calculation =================
for(h1 in 1:length(u)){
  for(i in 1:n.f){
    
    fpar = fpar.read[, -1]
    
    if(obs.dur[i] <= r[i]){ 
      
      #Caculate TaR forecast
      lam1 = fpar$lam1[i]
      del1 = fpar$del1[i]
      w1 = fpar$w1[i]
      b1 = besselK(w1, lam1/del1, expon.scaled = FALSE)
      b2 = besselK(w1, (lam1+1)/del1, expon.scaled = FALSE)
      
      tar = 2*u[h1]*(w1/2)^(lam1/del1) * b1
      z = findprob(myfunc1, interval=c(0,20), target=tar)
      
      fTaR[i] = fphi[i] * (2*z/w1)^(1/del1) * b1/b2
      
    } else {
      
      #Caculate TaR forecast
      lam2 = fpar$lam2[i]
      del2 = fpar$del2[i]
      w2 = fpar$w2[i]
      b1 = besselK(w2, lam2/del2, expon.scaled = FALSE)
      b2 = besselK(w2, (lam2+1)/del2, expon.scaled = FALSE)
      
      tar = 2*u[h1]*(w2/2)^(lam2/del2) * b1
      z = findprob(myfunc2, interval=c(0,20), target=tar)
      
      fTaR[i] = fphi[i] * (2*z/w2)^(1/del2) * b1/b2
      
    }
  }
  out.fTaR = cbind(out.fTaR, fTaR)
}

#Combine forecast value & forecast VaR
out.forc = data.frame(ObsDur=obs.dur[1:n.f], fphi=fphi[1:n.f], fTaR=out.fTaR)
colnames(out.forc) = c("ObsDur", "fphi", "fTaR_0.5", "fTaR_0.25", "fTaR_0.1", 
                       "fTaR_0.05", "fTaR_0.025", "fTaR_0.01")

write.csv(out.forc, "AAPL ObsDur, fphi & TaR Forecast_EGIG.TSCD.csv")


##############################################################################
############################################################################## 
# ================== Ratio of Violation Rate (RVR) ===================
I.VR = numeric(n.f)
VR = numeric(length(u))
RVR = numeric(length(u))

out.RVR = NULL

for(h2 in 1:length(u)){
  for(j in 1:n.f){
    
    I.VR[j] = ifelse(obs.dur[j] > out.fTaR[j, h2], 1, 0)
    
  }
  
  VR[h2] = mean(I.VR)
  RVR[h2] =  VR[h2] / (1-u[h2])
}

out.RVR = cbind(out.RVR, VR, RVR)
write.csv(out.RVR, "AAPL RVR Forecast_EGIG.TSCD.csv")


##############################################################################
##############################################################################
# ========== Backtesting: Kupiec, Dynamic Quantile, Christoffersen =============
out.fKup = NULL #Forecast Kupiec test
out.fKup.func = NULL #Forecast Kupiec buil-in func test
out.fDQ = NULL #Forecast Dynamic Quantile test
out.fCC = NULL #Forecast Christoffesen test

for(h3 in 1:length(u)){
  H0 = 'H0: Correct Exceedences'
  f.alpha = 0.05
  
  # ----------------- Kupiec Test -----------------
  f.Kup.p = 1-u[h3]
  f.exp.fail = n.f - round(u[h3]*n.f, 0) #expected failure
  
  f.fail = out.forc$ObsDur > out.forc[,(h3+2)] #failure
  
  f.num.fail = length(f.fail[f.fail==TRUE]) #num failure
  
  f.rat.fail = f.num.fail/n.f #x/T prop failure
  
  if(f.num.fail==0){
    f.last = 1
  } else {
    f.last = f.num.fail*log(f.rat.fail)
  }
  
  # Likelihood ratio statistic for Kupiec test
  f.LR = -2*(f.num.fail*log(f.Kup.p) + (n.f-f.num.fail)*log(1-f.Kup.p) - 
               f.last - (n.f-f.num.fail)*log(1-f.rat.fail))
  
  #critical value
  f.chiq = qchisq(p=1-f.alpha, df=1)
  
  #p-value
  f.p.value = 1-pchisq(q=f.LR, df=1)
  
  # Test decision
  if(f.p.value > f.alpha){
    f.con = 'Fail to Reject H0'
  }else{
    f.con = 'Reject H0'
  }
  
  out.fKup = rbind(out.fKup, data.frame(f.exp.fail=f.exp.fail,
                                        f.act.fail=f.num.fail,
                                        f.rat.fail=f.rat.fail,
                                        f.test.stat=f.LR,
                                        f.cri.v=f.chiq, f.pvalue=f.p.value,
                                        Decision=f.con))
  
  
  #Kupiec use build-in function ==============================================
  f.Kup.st.func = GAS::BacktestVaR(data=out.forc$ObsDur, VaR=out.forc[,(h3+2)], 
                                   alpha=u[h3])$LRuc[1]
  f.Kup.p.func = GAS::BacktestVaR(data=out.forc$ObsDur, VaR=out.forc[,(h3+2)], 
                                  alpha=u[h3])$LRuc[2]
  
  # Test decision
  if(f.Kup.p.func > f.alpha){
    f.conKup = 'Fail to Reject H0'
  }else{
    f.conKup = 'Reject H0'
  }
  out.fKup.func = rbind(out.fKup.func, data.frame(f.Kup.stat.func=f.Kup.st.func, 
                                                  f.Kup.pvalue.func=f.Kup.p.func, Dcision=f.conKup))
  #=============================================================================
  
  
  # ----------------- Dynamic Quantile Test -----------------
  f.DQ.st = GAS::BacktestVaR(data=out.forc$ObsDur, VaR=out.forc[,(h3+2)],
                             alpha=u[h3])$DQ$stat
  f.DQ.p = GAS::BacktestVaR(data=out.forc$ObsDur, VaR=out.forc[,(h3+2)],
                            alpha=u[h3])$DQ$pvalue
  
  if(f.DQ.p > f.alpha){
    f.conDQ = 'Fail to Reject H0'
  }else{
    f.conDQ = 'Reject H0'
  }
  out.fDQ = rbind(out.fDQ, data.frame(f.DQ.stat=f.DQ.st, f.DQ.pvalue=f.DQ.p,
                                      Dcision=f.conDQ))
  
  
  # ----------------- Christoffersen Test -----------------
  f.CC.st = GAS::BacktestVaR(data=out.forc$ObsDur, VaR=out.forc[,(h3+2)],
                             alpha=u[h3])$LRcc[1]
  f.CC.p = GAS::BacktestVaR(data=out.forc$ObsDur, VaR=out.forc[,(h3+2)],
                            alpha=u[h3])$LRcc[2]
  
  if(f.CC.p > f.alpha){
    f.conCC = 'Fail to Reject H0'
  }else{
    f.conCC = 'Reject H0'
  }
  out.fCC = rbind(out.fCC, data.frame(f.CC.stat=f.CC.st, f.CC.pvalue=f.CC.p,
                                      Dcision=f.conCC))
}


rownames(out.fKup) = c("up.risk_0.5", "up.risk_0.25", "up.risk_0.1",
                       "up.risk_0.05", "up.risk_0.025", "up.risk_0.01")

#Kupiec use build-in function ==============================================
rownames(out.fKup.func) = c("up.risk_0.5", "up.risk_0.25", "up.risk_0.1", 
                            "up.risk_0.05", "up.risk_0.025", "up.risk_0.01")

rownames(out.fDQ) = c("up.risk_0.5", "up.risk_0.25", "up.risk_0.1",
                      "up.risk_0.05", "up.risk_0.025", "up.risk_0.01")

rownames(out.fCC) = c("up.risk_0.5", "up.risk_0.25", "up.risk_0.1",
                      "up.risk_0.05", "up.risk_0.025", "up.risk_0.01")

write.csv(out.fKup, "AAPL Kupiec Forecast_EGIG.TSCD.csv")
write.csv(out.fKup.func, "AAPL Kupiec buil-in func Forecast_EGIG.TSCD.csv")
write.csv(out.fDQ, "AAPL DQ Forecast_EGIG.TSCD.csv")
write.csv(out.fCC, "AAPL CC Forecast_EGIG.TSCD.csv")
