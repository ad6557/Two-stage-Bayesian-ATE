library(BART)
PS.design = function(outcome,
                     treatment, #treatment vector
                     covardataset, #n x p covariate matrix where p is the number of covariates
                     bayes=TRUE, #draw from PS posterior distribution or estimate MLE
                     K=50){ #number of draws from the PS posterior distribution
  
  posteriors <- wbart(x.train = covardataset[treatment==1,], 
                      y.train = outcome[treatment==1],
                      x.test = covardataset, ndpost = K,nskip=5*K,printevery = 1000)
  if(bayes){
    score = t(posteriors[["yhat.test"]]) #f(x)
  } else {
    score = apply(posteriors[["yhat.test"]],2,mean)
  }
  
  return(score) #returns a n x K matrix
}