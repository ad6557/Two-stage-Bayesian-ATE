source("risk.data.generation.R")
source("risk.designfunction.R")
source("risk.analysisfunction1.R")

library(dplyr)

results=data.frame(matrix(ncol = 9, nrow = M))
colnames(results)=c("trueATE",
                    "BayesATE", "BayesSE", "BayesCI1", "BayesCI2", "B_coverage",
                    "BayesQ2.5", "BayesQ97.5", "Quantile_coverage")

model = "stan4bart"
for(m in 1:M){
  
  outcome=binarydata[[m]][["outcome"]]
  treatment=binarydata[[m]][["treatment"]]
  covardataset=binarydata[[m]][["covariates"]]
  trueATE=binarydata[[m]][["true"]]
  
  cat(m, "\n")

  K=100;S=1000; total=c()
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
  BayesATE = mean(stratATE,na.rm=TRUE)
  BayesSE = sqrt(mean(stratSE^2,na.rm=TRUE) + (1/K+1)*var(stratATE,na.rm=TRUE))
  BayesCI1 = BayesATE - 1.96*BayesSE
  BayesCI2 = BayesATE + 1.96*BayesSE
  
  B_coverage=between(trueATE,BayesCI1,BayesCI2)
  results[m,1:6]=data.frame(trueATE,BayesATE,BayesSE,BayesCI1,BayesCI2,B_coverage)
  
  Q1=quantile(total,prob=0.025,na.rm=TRUE);Q2=quantile(total,prob=0.975,na.rm=TRUE)
  results[m,7:9]=data.frame(Q1,Q2,between(trueATE,Q1,Q2))

}

results.final = na.omit(results)

mean(results.final$B_coverage)
mean(results.final$Quantile_coverage)
mean(results.final$trueATE-results.final$BayesATE)

save(beta0,beta1,D,binarydata,results,seeed,file=paste0(model,"_p",p,".RData"))
