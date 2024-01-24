library(BayesTree)
PS.design = function(outcome,
                     treatment, #treatment vector
                     covardataset, #n x p covariate matrix where p is the number of covariates
                     K){ #number of draws from the PS posterior distribution
  
  posteriors <- bart(x.train = covardataset[treatment==1,], y.train = outcome[treatment==1], 
                     x.test = covardataset, ndpost = K,verbose = FALSE)
  score = t(posteriors[["yhat.test"]]) #f(x)
  
  return(score) #returns a n x K matrix
}
