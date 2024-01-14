#library(BayesTree)
library(stan4bart)
PS.design = function(outcome,
                     treatment, #treatment vector
                     covardataset, #n x p covariate matrix where p is the number of covariates
                     K){ #number of draws from the PS posterior distribution
  
  #posteriors <- bart(x.train = covardataset[treatment==1,], y.train = outcome[treatment==1], 
  #                   x.test = covardataset, ndpost = K,verbose = FALSE)
  #score = t(posteriors[["yhat.test"]])
  
  dat = data.frame(outcome,treatment,covardataset)
  posteriors = stan4bart(outcome~bart(.-treatment)+treatment,dat,
                         chains = 1,seed = 0,iter = 2*K,warmup = K)
  score=extract(posteriors,type = "ev")
  
  return(score) #returns a n x K matrix
}
