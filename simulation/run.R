source("risk.data.generation.R")
source("risk.designfunction.R")
source("risk.analysisfunction.R")

rubinse = function(ates,sigmas){
  K = length(ates)
  sqrt(mean(sigmas^2,na.rm=TRUE) + var(ates,na.rm=TRUE))
}

bayes.output = function(stratATE,stratSE){
  BPSA.ATE = mean(stratATE,na.rm=TRUE)
  BPSA.lower = BPSA.ATE - 1.96*rubinse(stratATE,stratSE)
  BPSA.upper = BPSA.ATE + 1.96*rubinse(stratATE,stratSE)
  #cat("BPSA: ",round(BPSA.ATE,3), "(", round(BPSA.lower,3), ",", round(BPSA.upper,3), ")\n")
  return(c(BPSA.lower,BPSA.upper))
}

freq.output = function(PSA.ate,PSA.se){
  PSA.lower = PSA.ate - 1.96*PSA.se
  PSA.upper = PSA.ate + 1.96*PSA.se
  #cat("PSA: ",round(PSA.ate,3), "(", round(PSA.lower,3), ",", round(PSA.upper,3), ")\n")
  return(c(PSA.lower,PSA.upper))
}

library(dplyr)
coverage_BPSA=rep(NA,M)
coverage_PSA=rep(NA,M)

for(m in 1:M){
  
  outcome=continuousdata[[m]][["outcome"]]
  treatment=continuousdata[[m]][["treatment"]]
  covardataset=continuousdata[[m]][["covariates"]]
  trueATE=continuousdata[[m]][["benefitATE"]]
  
  #cat(m, "\n")
  #cat("true ATE: ", trueATE, "\n")
  ############ BPSA ##############
  K=800 # K=500 cannot converge
  BPSA.est.PS = PS.design(outcome,treatment,covardataset,bayes=TRUE,K=K)
  stratATE = rep(NA,K)
  stratSE = rep(NA,K)
  for(k in 1:K){
    strat = PS.analysis(BPSA.est.PS[,k],treatment,outcome,covardataset,bpsa=2,S=600) # S=400 cannot converge
    stratATE[k] = strat$ATE
    stratSE[k] = strat[["Average SE"]]
  }
  CI=bayes.output(stratATE,stratSE)
  coverage_BPSA[m]=between(trueATE,CI[1],CI[2])
  
  ############ PSA ##############
  PSA.est.PS = PS.design(outcome,treatment,covardataset,bayes=FALSE)
  strat = PS.analysis(PSA.est.PS,treatment,outcome,covardataset,bpsa = 1)
  
  PSA.ate = mean(strat$ATE)
  PSA.se = strat[["Average SE"]]
  
  CI=freq.output(PSA.ate,PSA.se)
  coverage_PSA[m]=between(trueATE,CI[1],CI[2])
  #cat("\n")
}

