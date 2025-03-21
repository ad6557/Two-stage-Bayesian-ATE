library(BART)
PS.design = function(outcome,
                     covardataset, #ncol:covariates, nrows:observations
                     id1, #training set
                     bayes=TRUE,
                     K=50){ #number of draws from the PBS posterior distribution
  
  posteriors <- wbart(x.train = covardataset[id1,], 
                      y.train = outcome[id1],
                      x.test = covardataset[-id1,], 
                      ndpost = K,nskip=5*K,printevery = 1000)
  
  if(bayes){
    score = t(posteriors[["yhat.test"]]) # n x K matrix
  } else {
    score = apply(posteriors[["yhat.test"]],2,mean) # vector of n
  }
  
  return(score)
}