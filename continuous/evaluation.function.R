library(beanz) # effect estimation for all subgroups in one function

PS.analysis = function(psvector, #vector of estimated PBS
                       treatment,  #treatment vector
                       outcome,   #outcome vector 
                       covardataset,  #ncol:covariates, nrows:observations
                       minCond){
  #data structure creation
  dataset = data.frame(outcome,treatment,psvector)
  dataset$indicator = ifelse(dataset$psvector>minCond,2,1)
  
  subgrp.effect <- bzGetSubgrpRaw(dataset,resptype = "continuous", 
                                  var.resp = "outcome",var.trt = "treatment",
                                  var.cov = "indicator",var.censor = "psvector")

  or.all <- glm(outcome~treatment, family=gaussian);
  eff    <- coef(summary(or.all))[2, 1];
  eff.sd <- coef(summary(or.all))[2, 2];
  
  estimated.effect = c(subgrp.effect$Estimate,eff)
  estimated.se = c(sqrt(subgrp.effect$Variance),eff.sd)
  return(list("ATE" = estimated.effect,"SE" = estimated.se))
}
