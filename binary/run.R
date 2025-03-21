source("data.generation.R")
source("design.function.R")
source("evaluation.function.R")

library(dplyr)

################################## functions ##################################
rubin.se = function(stratATE,stratSE,K){
  BayesSE = sqrt(mean(stratSE^2,na.rm=TRUE) + (1/K+1)*var(stratATE,na.rm=TRUE))
  return(BayesSE)
}

output = function(ATE,SE,trueATE){
  CI1 = ATE - 1.96*SE
  CI2 = ATE + 1.96*SE
  C_coverage=between(trueATE,CI1,CI2)
  return(data.frame(SE,CI1,CI2,C_coverage))
}

################################## run model ##################################
results_NLR=data.frame(matrix(ncol = 10, nrow = M))
colnames(results_NLR)=c("trueATE","EstATE", 
                        "NaiveSE", "NaiveCI1", "NaiveCI2", "F_coverage",
                        "BayesSE", "BayesCI1", "BayesCI2", "B_coverage")
results_LR = results_NLR
results_all = results_NLR
K=100

results_NLR[,1]=log_NLR
results_LR[,1]=log_LR
results_all[,1]=log_all

for(m in 1:M){
  
  outcome=binarydata[[m]][["outcome"]]
  treatment=binarydata[[m]][["treatment"]]
  covardataset=binarydata[[m]][["covariates"]]
  id1 = sample((n/2+1):n,n/4)
  
  cat("m",m, "\n")
  
  # Naive model
  PSA.est.PS = PS.design(outcome,covardataset,id1,bayes=FALSE)
  MLE = PS.analysis(PSA.est.PS,treatment[-id1],outcome[-id1],covardataset[-id1,],minCond)
  # point estimate
  results_NLR[m,2]=MLE$ATE[1]
  results_LR[m,2]=MLE$ATE[2]
  results_all[m,2]=MLE$ATE[3]
  # CI and coverage
  results_NLR[m,3:6]=output(results_NLR[m,2],MLE$SE[1],log_NLR)
  results_LR[m,3:6]=output(results_LR[m,2],MLE$SE[2],log_LR)
  results_all[m,3:6]=output(results_all[m,2],MLE$SE[3],log_all)
  
  # Two Stage model
  total_NLR=c();total_LR=c();total_all=c();
  BPSA.est.PS = PS.design(outcome,covardataset,id1,bayes=TRUE,K=K)
  MLE_ATE = matrix(NA,K,3) # NLR LR ALL
  MLE_SE = matrix(NA,K,3) # NLR LR ALL
  for(k in 1:K){
    psvector = BPSA.est.PS[,k]
    MLE = PS.analysis(psvector,treatment[-id1],outcome[-id1],covardataset[-id1,],minCond)
    MLE_ATE[k,] = MLE$ATE
    MLE_SE[k,] = MLE$SE
  }
  # uncertainty
  se_NLR = rubin.se(MLE_ATE[,1],MLE_SE[,1],K)
  se_LR = rubin.se(MLE_ATE[,2],MLE_SE[,2],K)
  se_all = rubin.se(MLE_ATE[,3],MLE_SE[,3],K)
  # CI and coverage
  results_NLR[m,7:10]=output(results_NLR[m,2],se_NLR,log_NLR)
  results_LR[m,7:10]=output(results_LR[m,2],se_LR,log_LR)
  results_all[m,7:10]=output(results_all[m,2],se_all,log_all)
}


################################## results ##################################

# NLR
bias_F_NLR = mean(results_NLR$trueATE-results_NLR$EstATE)
bias_B_NLR = mean(results_NLR$trueATE-results_NLR$EstATE)

var_F_NLR = var(results_NLR$EstATE)
var_B_NLR = var(results_NLR$EstATE)

MSE_F_NLR = bias_F_NLR^2+var_F_NLR
MSE_B_NLR = bias_B_NLR^2+var_B_NLR

se_F_NLR = mean(results_NLR$NaiveSE)
se_B_NLR = mean(results_NLR$BayesSE)

cvrg_F_NLR = mean(results_NLR$F_coverage)
cvrg_B_NLR = mean(results_NLR$B_coverage)

# LR
bias_F_LR = mean(results_LR$trueATE-results_LR$EstATE)
bias_B_LR = mean(results_LR$trueATE-results_LR$EstATE)

var_F_LR = var(results_LR$EstATE)
var_B_LR = var(results_LR$EstATE)

MSE_F_LR = bias_F_LR^2+var_F_LR
MSE_B_LR = bias_B_LR^2+var_B_LR

se_F_LR = mean(results_LR$NaiveSE)
se_B_LR = mean(results_LR$BayesSE)

cvrg_F_LR = mean(results_LR$F_coverage)
cvrg_B_LR = mean(results_LR$B_coverage)

# all
bias_F_all = mean(results_all$trueATE-results_all$EstATE)
bias_B_all = mean(results_all$trueATE-results_all$EstATE)

var_F_all = var(results_all$EstATE)
var_B_all = var(results_all$EstATE)

MSE_F_all = bias_F_all^2+var_F_all
MSE_B_all = bias_B_all^2+var_B_all

se_F_all = mean(results_all$NaiveSE)
se_B_all = mean(results_all$BayesSE)

cvrg_F_all = mean(results_all$F_coverage)
cvrg_B_all = mean(results_all$B_coverage)

# arrange output table
overall_table <- data.frame(MSE__F = c(MSE_F_NLR,MSE_F_LR,MSE_F_all),
                            var_F = c(var_F_NLR,var_F_LR,var_F_all),
                            Avg_Bias_F=c(bias_F_NLR,bias_F_LR,bias_F_all),
                            Avg_SE_F=c(se_F_NLR,se_F_LR,se_F_all),
                            coverage_F=c(cvrg_F_NLR,cvrg_F_LR,cvrg_F_all),
                            
                            MSE__B = c(MSE_B_NLR,MSE_B_LR,MSE_B_all),
                            var_B = c(var_B_NLR,var_B_LR,var_B_all),
                            Avg_Bias_B=c(bias_B_NLR,bias_B_LR,bias_B_all),
                            Avg_SE_B=c(se_B_NLR,se_B_LR,se_B_all),
                            coverage_B=c(cvrg_B_NLR,cvrg_B_LR,cvrg_B_all),
                            
                            row.names=c('NLR', 'LR', 'all'))
# print
apply(overall_table, 2, function(x) round(x, 3))

