#required packages 
library(stan4bart)
#analysis function
PS.analysis = function(psvector, ##psvector = input vector of ps 
                       treatment,  ##treatment = input indicator vector of treated/control
                       outcome,   ##outcome = input vector of outcomes
                       covardataset,  ##covardataset = matrix of covariates, ncol = # of covariates, nrows = # of observations
                       model, # bartc or stan4bart
                       S = 500){ ##draws from the posterior distribution of delta
  #data structure creation
  dat = data.frame(outcome,treatment,covardataset)
  cutoff = 0.5 # if stage1 uses bart() then cutoff at 0
  dat$indicator = ifelse(psvector>cutoff,0,1) 
  
  if (model=="stan4bart") {
    
    fit <- stan4bart(outcome ~ bart(.-indicator)+indicator, dat,
                     treatment = treatment,chains = 1,seed = 0,iter = 2*S,warmup = S)
    ev = extract(fit,type="ev")
    ev0 = ev[dat$indicator==0,]; treatment0 = treatment[dat$indicator==0];
    p1.0 <- apply(ev0[which(treatment0==1),],2,mean)
    p0.0 <- apply(ev0[which(treatment0==0),],2,mean)
    ates = p1.0-p0.0
    
  }
  
  estimated.effect = mean(ates)
  estimated.se = sd(ates)
       
  return(list("ATE" = estimated.effect,"Average SE" = estimated.se,"ATES" = as.vector(ates)))
}
