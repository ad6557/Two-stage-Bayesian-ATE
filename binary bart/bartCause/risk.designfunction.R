
library(bartCause)

PS.design = function(outcome,
                     treatment, #treatment vector
                     covardataset, #n x p covariate matrix where p is the number of covariates
                     K){ #number of draws from the PS posterior distribution
  
  posteriors <- bartc(outcome, treatment, covardataset, 
                      method.trt = "none", estimand = "ate",
                      n.samples=K, n.burn = 10*K, n.chains = 1,verbose=FALSE)
  score = extract(posteriors,type = "mu.1")
  return(t(score)) #returns a n x K matrix
}
