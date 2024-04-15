

model="linear"
M = 200
n = 200
p = 10
minCond = 0
source("risk.designfunction.R")
source("risk.analysisfunction.R")

library(dplyr)

################################## functions ##################################
bayes.output = function(stratATE,stratSE,K,trueATE){
  BayesATE = mean(stratATE,na.rm=TRUE)
  BayesSE = sqrt(mean(stratSE^2,na.rm=TRUE) + (1/K+1)*var(stratATE,na.rm=TRUE))
  BayesCI1 = BayesATE - 1.96*BayesSE
  BayesCI2 = BayesATE + 1.96*BayesSE
  B_coverage=between(trueATE,BayesCI1,BayesCI2)
  return(data.frame(BayesATE,BayesSE,BayesCI1,BayesCI2,B_coverage))
}

freq.output = function(FreqATE,FreqSE,trueATE){
  FreqCI1 = FreqATE - 1.96*FreqSE
  FreqCI2 = FreqATE + 1.96*FreqSE
  F_coverage=between(trueATE,FreqCI1,FreqCI2)
  return(data.frame(trueATE,FreqATE,FreqSE,FreqCI1,FreqCI2,F_coverage))
}

bayes.quantile.output = function(total,trueATE){
  Q1=quantile(total,prob=0.025,na.rm=TRUE)
  Q2=quantile(total,prob=0.975,na.rm=TRUE)
  Q_coverage=between(trueATE,Q1,Q2)
  return(data.frame(mean(total),sd(total),Q1,Q2,Q_coverage))
}

new_NLR=data.frame(matrix(ncol = 21, nrow = 200))
colnames(new_NLR)=c("trueATE",
                        "FreqATE", "FreqSE", "FreqCI1", "FreqCI2", "F_coverage",
                        "BayesATE", "BayesSE", "BayesCI1", "BayesCI2", "B_coverage",
                        "BayesQATE","BayesQSE","BayesQ2.5", "BayesQ97.5", "Quantile_coverage",
                        "BayesRATE","BayesRSE","BayesRCI1", "BayesRCI2", "B_R_coverage")
new_LR = new_NLR
new_all = new_NLR

minCond = 0
for(m in 1:200){
  
  outcome=binarydata[[m]][["outcome"]]
  treatment=binarydata[[m]][["treatment"]]
  covardataset=binarydata[[m]][["covariates"]]
  
  cat("m",m, "\n")
  
  # bayes
  K=100; 
  total_NLR=c();total_LR=c();total_all=c();
  BPSA.est.PS = PS.design(outcome,treatment,covardataset,bayes=TRUE,K=K)
  MLE_ATE = matrix(NA,K,3); BAY_ATE = MLE_ATE # NLR LR ALL
  MLE_SE = matrix(NA,K,3); BAY_SE = MLE_SE # NLR LR ALL
  for(k in 1:K){
    psvector = BPSA.est.PS[,k]
    MLE = PS.analysis(psvector,treatment,outcome,covardataset,minCond,bayes=FALSE)
    MLE_ATE[k,] = MLE$ATE
    MLE_SE[k,] = MLE$SE
  }
  
  new_NLR[m,7:11]=bayes.output(MLE_ATE[,1],MLE_SE[,1],K,log(OR_NLR))
  new_LR[m,7:11]=bayes.output(MLE_ATE[,2],MLE_SE[,2],K,log(OR_LR))
  new_all[m,7:11]=bayes.output(MLE_ATE[,3],MLE_SE[,3],K,log(OR_all))
  
  # freq
  PSA.est.PS = PS.design(outcome,treatment,covardataset,bayes=FALSE)
  MLE = PS.analysis(PSA.est.PS,treatment,outcome,covardataset,minCond,bayes=FALSE)
  
  new_NLR[m,1:6]=freq.output(MLE$ATE[1],MLE$SE[1],log(OR_NLR))
  new_LR[m,1:6]=freq.output(MLE$ATE[2],MLE$SE[2],log(OR_LR))
  new_all[m,1:6]=freq.output(MLE$ATE[3],MLE$SE[3],log(OR_all))
}

results_NLR = results_NLR[1:200,]
results_LR = results_LR[1:200,]
results_all = results_all[1:200,]

results_NLR[,1:11] = new_NLR[,1:11]
results_LR[,1:11] = new_LR[,1:11]
results_all[,1:11] = new_all[,1:11]


# save
if (dim(binarydata[[1]][["covariates"]])[2] == p){
  save(results_NLR,results_LR,results_all,binarydata,OR_NLR,OR_LR,OR_all,
       file=paste0(model,"_n",n,"_p",p,"_M",M,"minCond",minCond,".RData"))
} else {
  save(results_NLR,results_LR,results_all,binarydata,OR_NLR,OR_LR,OR_all,
       file=paste0(model,"_n",n,"_p",p,"actually5","_M",M,"minCond",minCond,".RData"))
}

