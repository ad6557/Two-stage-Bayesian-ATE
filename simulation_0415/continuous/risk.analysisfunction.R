#required packages 
require(rstanarm)
library(beanz)
#analysis function
PS.analysis = function(psvector, ##psvector = input vector of ps 
                       treatment,  ##treatment = input indicator vector of treated/control
                       outcome,   ##outcome = input vector of outcomes
                       covardataset,  ##ncol = #covariates, nrows = #observations
                       minCond,
                       bayes=TRUE){
  #data structure creation
  dataset = data.frame(outcome,treatment,psvector)
  dataset$indicator = ifelse(dataset$psvector>minCond,2,1)
  
  ############# subgroup effect estimation #############
  subgrp.effect <- bzGetSubgrpRaw(dataset,resptype = "continuous", 
                                  var.resp = "outcome",var.trt = "treatment",
                                  var.cov = "indicator",var.censor = "psvector")
  # "var.censor" doesn't influence results
  or.all <- glm(outcome~treatment, family=gaussian);
  eff    <- coef(summary(or.all))[2, 1];
  eff.sd <- coef(summary(or.all))[2, 2];
  
  if(bayes==FALSE){
    
    estimated.effect = c(subgrp.effect$Estimate,eff)
    estimated.se = c(sqrt(subgrp.effect$Variance),eff.sd) # var not sd, need sqrt
    return(list("ATE" = estimated.effect,"SE" = estimated.se))
    
  } else {
    
    var.estvar <- c("Estimate", "Variance")
    
    suppressWarnings(
      rst.fs <- bzCallStan("fs", dat.sub=subgrp.effect,
                           var.estvar = var.estvar, var.cov = "indicator",
                           par.pri = c(B=10000), chains=2, iter=2000,
                           warmup=1000, seed=1000, cores=1)
    )
    suppressWarnings(
      rst.nse <- bzCallStan("nse", dat.sub=subgrp.effect,
                            var.estvar = var.estvar, var.cov = "indicator",
                            par.pri = c(B=10000), chains=2, iter=2000,
                            warmup=1000, seed=1000, cores=1)
    )
    tbl.sub <- bzSummary(rst.fs, ref.stan.rst=rst.nse)
    mean_NLR = as.numeric(tbl.sub[1,2])
    mean_LR = as.numeric(tbl.sub[2,2])
    mean_all = as.numeric(tbl.sub[3,2])
    sd_NLR = as.numeric(tbl.sub[1,3])
    sd_LR = as.numeric(tbl.sub[2,3])
    sd_all = as.numeric(tbl.sub[3,3])
    
    ates_NLR = rst.fs[["smps"]][,,1]
    ates_LR = rst.fs[["smps"]][,,2]
    ates_all = rst.nse[["smps"]][,,1]
    return(list("ATES" = list(ates_NLR,ates_LR,ates_all),
                "MEAN" = c(mean_NLR,mean_LR,mean_all),
                "SD" = c(sd_NLR,sd_LR,sd_all)))
  }
  
}
