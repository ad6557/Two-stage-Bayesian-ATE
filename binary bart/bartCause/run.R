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

results=data.frame(matrix(ncol = 9, nrow = M))
colnames(results)=c("trueATE",
                    "BayesATE", "BayesSE", "BayesCI1", "BayesCI2", "BPSAcoverage",
                    "BayesQ2.5", "BayesQ97.5", "BPSAQcoverage")

model = "bartc_interaction"
for(m in 1:M){
  
  outcome=binarydata[[m]][["outcome"]]
  treatment=binarydata[[m]][["treatment"]]
  covardataset=binarydata[[m]][["covariates"]]
  trueATE=binarydata[[m]][["true"]]
  
  cat(m, "\n")

  K=50;S=1000; total=c()
  BPSA.est.PS = PS.design(outcome,treatment,covardataset,K=K)
  stratATE = rep(NA,K)
  stratSE = rep(NA,K)
  for(k in 1:K){
    psvector = BPSA.est.PS[,k]
    strat = PS.analysis(psvector,treatment,outcome,covardataset,model,S=S)
    stratATE[k] = strat$ATE
    stratSE[k] = strat[["Average SE"]]
    total = c(total,strat$ATES)
  }
  CI=bayes.output(stratATE,stratSE,K)
  coverage_BPSA=between(trueATE,CI[1],CI[2])
  results[m,1:6]=data.frame(trueATE,CI[3],CI[4],CI[1],CI[2],coverage_BPSA)
  
  Q1=quantile(total,prob=0.025,na.rm=TRUE);Q2=quantile(total,prob=0.975,na.rm=TRUE)
  results[m,7:9]=data.frame(Q1,Q2,between(trueATE,Q1,Q2))

}

results.final = na.omit(results)

mean(results.final$BPSAQcoverage)
mean(results.final$BPSAcoverage)
mean(results.final$trueATE-results.final$BayesATE)

save(beta0,beta1,D,binarydata,results,seeed,file=paste0(model,"_p",p,"bart.RData"))
