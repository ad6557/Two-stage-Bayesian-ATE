source("risk.data.generation.R")
source("risk.designfunction.R")
source("risk.analysisfunction1.R")

bayes.output = function(stratATE,stratSE,K){
  BPSA.ATE = mean(stratATE,na.rm=TRUE)
  BPSA.SE = sqrt(mean(stratSE^2,na.rm=TRUE) + (1/K+1)*var(stratATE,na.rm=TRUE))
  BPSA.lower = BPSA.ATE - 1.96*BPSA.SE
  BPSA.upper = BPSA.ATE + 1.96*BPSA.SE
  #cat("BPSA: ",round(BPSA.ATE,3), "(", round(BPSA.lower,3), ",", round(BPSA.upper,3), ")\n")
  return(c(BPSA.lower,BPSA.upper,BPSA.ATE,BPSA.SE))
}

freq.output = function(PSA.ate,PSA.se){
  PSA.lower = PSA.ate - 1.96*PSA.se
  PSA.upper = PSA.ate + 1.96*PSA.se
  #cat("PSA: ",round(PSA.ate,3), "(", round(PSA.lower,3), ",", round(PSA.upper,3), ")\n")
  return(c(PSA.lower,PSA.upper))
}

library(dplyr)
#coverage_BPSA=rep(NA,M)
#coverage_PSA=rep(NA,M)
results=data.frame(matrix(ncol = 14, nrow = M))
colnames(results)=c("trueATE",
                    "BayesATE", "BayesSE", "BayesCI1", "BayesCI2", "BPSAcoverage",
                    "FreqATE", "FreqSE", "FreqCI1", "FreqCI2", "PSAcoverage", 
                    "BayesQ2.5", "BayesQ97.5", "BPSAQcoverage")

for(m in 1:M){
  
  outcome=continuousdata[[m]][["outcome"]]
  treatment=continuousdata[[m]][["treatment"]]
  covardataset=continuousdata[[m]][["covariates"]]
  trueATE=continuousdata[[m]][["benefitATE"]]
  
  cat(m, "\n")
  #cat("true ATE: ", trueATE, "\n")
  ############ BPSA ##############
  K=50;S=1000; total=c()
  BPSA.est.PS = PS.design(outcome,treatment,covardataset,bayes=TRUE,K=K)
  stratATE = rep(NA,K)
  stratSE = rep(NA,K)
  for(k in 1:K){
    strat = PS.analysis(BPSA.est.PS[,k],treatment,outcome,covardataset,bpsa=2,S=S)
    stratATE[k] = strat$ATE
    stratSE[k] = strat[["Average SE"]]
    total = c(total,strat$ATES)
  }
  CI=bayes.output(stratATE,stratSE,K)
  coverage_BPSA=between(trueATE,CI[1],CI[2])
  results[m,1:6]=data.frame(trueATE,CI[3],CI[4],CI[1],CI[2],coverage_BPSA)
  
  Q1=quantile(total,prob=0.025);Q2=quantile(total,prob=0.975)
  results[m,12:14]=data.frame(Q1,Q2,between(trueATE,Q1,Q2))
  ############ PSA ##############
  PSA.est.PS = PS.design(outcome,treatment,covardataset,bayes=FALSE)
  strat = PS.analysis(PSA.est.PS,treatment,outcome,covardataset,bpsa = 1)
  
  PSA.ate = mean(strat$ATE)
  PSA.se = strat[["Average SE"]]
  
  CI=freq.output(PSA.ate,PSA.se)
  coverage_PSA=between(trueATE,CI[1],CI[2])
  results[m,7:11]=data.frame(PSA.ate,PSA.se,CI[1],CI[2],coverage_PSA)
  #cat("\n")
}

mean(results$BPSAQcoverage)
mean(results$BPSAcoverage)
mean(results$PSAcoverage)
mean(results$trueATE-results$BayesATE)
mean(results$trueATE-results$FreqATE)

save(list(results,beta0,beta1),file=paste0("seed",seeed,"M",M,"p",p,"_interaction.RData"))
