##Required functions and packages
expit = function(x){
  exp(x)/(1+exp(x))
}

require(rstanarm) #require(MCMCpack)

###Design function 
PS.design = function(outcome,
                     treatment, #treatment vector
                     covardataset, #n x p covariate matrix where p is the number of covariates
                     bayes=TRUE, #draw from PS posterior distribution or estimate MLE
                     K){ #number of draws from the PS posterior distribution
  
  dataset = data.frame(outcome=outcome,covardataset=covardataset,treatment=treatment)
  n = dim(covardataset)[1] #sample size
  p = dim(covardataset)[2]
  thin = 10
  
  ###Bayesian PS estimation
  if(bayes){
    #posteriors = as.matrix(MCMCregress(outcome~treatment*covardataset,data=dataset,
    #                                   burnin=2*K,mcmc=K*100,thin=100))
    # beta1 = posteriors[,c(2,(2+p+1):(2+p+5))]
    posteriors <- stan_glm(outcome~treatment*covardataset,data=dataset,
                           family = gaussian(link = "identity"),refresh=0,
                           iter = (thin+2)*K,warmup = 2*K,thin = 4*thin) # 4 chains
    beta1 = as.matrix(posteriors)[,c(2,(2+p+1):(2+p+5))]
    covars = as.matrix(cbind(rep(1,n),covardataset))
    score = covars%*%t(beta1) # benefit indicator = [1 X]*beta1
  } else {
    #########Frequentist PS estimation
    fit = lm(outcome~treatment*covardataset,data=dataset)
    beta1 = as.matrix(fit$coefficients)[c(2,8:12),]
    covars = as.matrix(cbind(rep(1,n),covardataset))
    score = covars%*%beta1 # benefit indicator = [1 X]*beta1
  }
  
    return(score) #returns a n x K matrix for BPSA, n x 1 for frequentist PSA
}
