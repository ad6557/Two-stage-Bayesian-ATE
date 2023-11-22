#required packages 
require(rstanarm)

#analysis function
PS.analysis = function(psvector, ##psvector = input vector of ps 
                       treatment,  ##treatment = input indicator vector of treated/control
                       outcome,   ##outcome = input vector of outcomes
                       covardataset,  ##covardataset = matrix of covariates, ncol = # of covariates, nrows = # of observations
                       bpsa = 1, ##bpsa-1 uses asymptotic methods (+ sandwich estimator), bpsa-2 samples from the conditional posterior distribution of delta
                       S = 500){ ##draws from the posterior distribution of delta
  #data structure creation
  dataset = data.frame(outcome,treatment,psvector,covardataset)
  dataset = dataset[dataset$psvector>0,]
  
  if(bpsa==1){# BPSA-1 : MLE linear regression with interactions
    fit = lm(outcome~treatment,data=dataset) # outcome~.-psvector / outcome~treatment
    estimated.effect = summary(fit)$coef[2,1]
    estimated.se = summary(fit)$coef[2,2]
    
    } else{#perform bayesian linear regression 
      tryCatch({
        #posteriors = as.data.frame(MCMCregress(outcome~.-psvector,data=dataset,
        #                                       burnin=S*5,mcmc=S*50,thin=50))
        posteriors <- stan_glm(outcome~treatment,data=dataset, # outcome~.-psvector / outcome~treatment
                               family = gaussian(link = "identity"),refresh=0,prior = NULL)
        draws = sample(c(1:4000),S) 
        ates = as.matrix(posteriors)[draws,"treatment"]
        estimated.effect = mean(ates)
        estimated.se = sd(ates)

      } , error=function(e){})
      
    }

  return(list("ATE" = estimated.effect,"Average SE" = estimated.se))
}
