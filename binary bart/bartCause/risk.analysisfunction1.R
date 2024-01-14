
#analysis function
PS.analysis = function(psvector, ##psvector = input vector of ps 
                       treatment,  ##treatment = input indicator vector of treated/control
                       outcome,   ##outcome = input vector of outcomes
                       covardataset,  ##covardataset = matrix of covariates, ncol = # of covariates, nrows = # of observations
                       model, # bartc or stan4bart
                       S = 500){ ##draws from the posterior distribution of delta
  #data structure creation
  dat = data.frame(outcome,treatment,covardataset)
  dat$indicator = ifelse(psvector>0.5,0,1)
  
  if(model=="bartc_interaction"){
    
    fit <- bartc(outcome, treatment, covardataset, group.by = dat$indicator, group.effects = TRUE,
                 n.samples = S, n.burn = S, n.chains = 2,verbose = FALSE,seed = 0)
    summary(fit)
    cate = extract(fit,type = "cate")
    ates = cate[["0"]]
    
  } else if (model=="bartc_subgroup") {
    
    fit <- bartc(outcome[dat$indicator==0], treatment[dat$indicator==0], covardataset[dat$indicator==0,],
                 n.samples = S, n.burn = S, n.chains = 1,verbose = FALSE)
    summary(fit) 
    ates = extract(fit,type = "cate")
    
  }
  
  estimated.effect = mean(ates)
  estimated.se = sd(ates)
       
  return(list("ATE" = estimated.effect,"Average SE" = estimated.se,"ATES" = as.vector(ates)))
}
