model="nonlinear"

source(paste0(model,".data.generation.R"))
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

################################## run model ##################################
results_NLR=data.frame(matrix(ncol = 21, nrow = M))
colnames(results_NLR)=c("trueATE",
                        "FreqATE", "FreqSE", "FreqCI1", "FreqCI2", "F_coverage",
                        "BayesATE", "BayesSE", "BayesCI1", "BayesCI2", "B_coverage",
                        "BayesQATE","BayesQSE","BayesQ2.5", "BayesQ97.5", "Quantile_coverage",
                        "BayesRATE","BayesRSE","BayesRCI1", "BayesRCI2", "B_R_coverage")
results_LR = results_NLR
results_all = results_NLR

for(m in 1:M){
  
  outcome=continuousdata[[m]][["outcome"]]
  treatment=continuousdata[[m]][["treatment"]]
  covardataset=continuousdata[[m]][["covariates"]]
  
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
    BAYE= PS.analysis(psvector,treatment,outcome,covardataset,minCond,bayes=TRUE)
    total_NLR = c(total_NLR,BAYE[["ATES"]][[1]])
    total_LR = c(total_LR,BAYE[["ATES"]][[2]])
    total_all = c(total_all,BAYE[["ATES"]][[3]])
    BAY_ATE[k,] = BAYE$MEAN
    BAY_SE[k,] = BAYE$SD
  }
  
  results_NLR[m,7:11]=bayes.output(MLE_ATE[,1],MLE_SE[,1],K,ATE_NLR)
  results_LR[m,7:11]=bayes.output(MLE_ATE[,2],MLE_SE[,2],K,ATE_LR)
  results_all[m,7:11]=bayes.output(MLE_ATE[,3],MLE_SE[,3],K,ATE_all)
  
  results_NLR[m,12:16]=bayes.quantile.output(total_NLR,ATE_NLR)
  results_LR[m,12:16]=bayes.quantile.output(total_LR,ATE_LR)
  results_all[m,12:16]=bayes.quantile.output(total_all,ATE_all)
  
  results_NLR[m,17:21]=bayes.output(BAY_ATE[,1],BAY_SE[,1],K,ATE_NLR)
  results_LR[m,17:21]=bayes.output(BAY_ATE[,2],BAY_SE[,2],K,ATE_LR)
  results_all[m,17:21]=bayes.output(BAY_ATE[,3],BAY_SE[,3],K,ATE_all)
  
  # freq
  PSA.est.PS = PS.design(outcome,treatment,covardataset,bayes=FALSE)
  MLE = PS.analysis(PSA.est.PS,treatment,outcome,covardataset,minCond,bayes=FALSE)
  
  results_NLR[m,1:6]=freq.output(MLE$ATE[1],MLE$SE[1],ATE_NLR)
  results_LR[m,1:6]=freq.output(MLE$ATE[2],MLE$SE[2],ATE_LR)
  results_all[m,1:6]=freq.output(MLE$ATE[3],MLE$SE[3],ATE_all)
}

# save
if (dim(continuousdata[[1]][["covariates"]])[2] == p){
  save(results_NLR,results_LR,results_all,continuousdata,ATE_NLR,ATE_LR,ATE_all,
       file=paste0(model,"_n",n,"_p",p,"_M",M,"minCond",minCond,".RData"))
} else {
  save(results_NLR,results_LR,results_all,continuousdata,ATE_NLR,ATE_LR,ATE_all,
       file=paste0(model,"_n",n,"_p",p,"actually5","_M",M,"minCond",minCond,".RData"))
}

################################## results ##################################

ind_LR = unique(c(which(results_LR$FreqSE>10), which(results_LR$BayesSE>10)))
ind_NLR = unique(c(which(results_NLR$FreqSE>10), which(results_NLR$BayesSE>10)))

#results_LR = results_LR[-ind_LR,]
#results_NLR = results_NLR[-ind_NLR,]

# NLR
bias_F_NLR = mean(results_NLR$trueATE-results_NLR$FreqATE)
bias_B_NLR = mean(results_NLR$trueATE-results_NLR$BayesATE)
bias_Q_NLR = mean(results_NLR$trueATE-results_NLR$BayesQATE)
bias_BR_NLR = mean(results_NLR$trueATE-results_NLR$BayesRATE)

se_F_NLR = mean(results_NLR$FreqSE)
se_B_NLR = mean(results_NLR$BayesSE)
se_Q_NLR = mean(results_NLR$BayesQSE)
se_BR_NLR = mean(results_NLR$BayesRSE)

cvrg_F_NLR = mean(results_NLR$F_coverage)
cvrg_B_NLR = mean(results_NLR$B_coverage)
cvrg_Q_NLR = mean(results_NLR$Quantile_coverage)
cvrg_BR_NLR = mean(results_NLR$B_R_coverage)

# LR
bias_F_LR = mean(results_LR$trueATE-results_LR$FreqATE)
bias_B_LR = mean(results_LR$trueATE-results_LR$BayesATE)
bias_Q_LR = mean(results_LR$trueATE-results_LR$BayesQATE)
bias_BR_LR = mean(results_LR$trueATE-results_LR$BayesRATE)

se_F_LR = mean(results_LR$FreqSE)
se_B_LR = mean(results_LR$BayesSE)
se_Q_LR = mean(results_LR$BayesQSE)
se_BR_LR = mean(results_LR$BayesRSE)

cvrg_F_LR = mean(results_LR$F_coverage)
cvrg_B_LR = mean(results_LR$B_coverage)
cvrg_Q_LR = mean(results_LR$Quantile_coverage)
cvrg_BR_LR = mean(results_LR$B_R_coverage)

# all
bias_F_all = mean(results_all$trueATE-results_all$FreqATE)
bias_B_all = mean(results_all$trueATE-results_all$BayesATE)
bias_Q_all = mean(results_all$trueATE-results_all$BayesQATE)
bias_BR_all = mean(results_all$trueATE-results_all$BayesRATE)

se_F_all = mean(results_all$FreqSE)
se_B_all = mean(results_all$BayesSE)
se_Q_all = mean(results_all$BayesQSE)
se_BR_all = mean(results_all$BayesRSE)

cvrg_F_all = mean(results_all$F_coverage)
cvrg_B_all = mean(results_all$B_coverage)
cvrg_Q_all = mean(results_all$Quantile_coverage)
cvrg_BR_all = mean(results_all$B_R_coverage)

# arrange output table
overall_table <- data.frame(Avg_Bias_F=c(bias_F_NLR,bias_F_LR,bias_F_all),
                            Avg_SE_F=c(se_F_NLR,se_F_LR,se_F_all),
                            coverage_F=c(cvrg_F_NLR,cvrg_F_LR,cvrg_F_all),
                            
                            Avg_Bias_B=c(bias_B_NLR,bias_B_LR,bias_B_all),
                            Avg_SE_B=c(se_B_NLR,se_B_LR,se_B_all),
                            coverage_B=c(cvrg_B_NLR,cvrg_B_LR,cvrg_B_all),
                            
                            Avg_Bias_Q=c(bias_Q_NLR,bias_Q_LR,bias_Q_all),
                            Avg_SE_Q=c(se_Q_NLR,se_Q_LR,se_Q_all),
                            coverage_Q=c(cvrg_Q_NLR,cvrg_Q_LR,cvrg_Q_all),
                            
                            Avg_Bias_R=c(bias_BR_NLR,bias_BR_LR,bias_BR_all),
                            Avg_SE_R=c(se_BR_NLR,se_BR_LR,se_BR_all),
                            coverage_R=c(cvrg_BR_NLR,cvrg_BR_LR,cvrg_BR_all),
                            
                            row.names=c('NLR', 'LR', 'all'))
# print
apply(overall_table, 2, function(x) round(x, 3))


