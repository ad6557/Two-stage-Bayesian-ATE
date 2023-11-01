#load packages + required functions
require(MASS)

expit = function(x){
  exp(x)/(1 + exp(x))
}

###########
n =100 #sample size
M = 500 #replications

p = 5
beta0 = matrix(c(0.1,seq(0.2,1,length=p)),nrow=1) #covariates for PS model
beta1 = matrix(c(0.3,c(-0.8,-0.5,-0.1,0.5,1)),nrow=1)
D = diag(rep(1,p)) #correlation matrix for covaritates: all covariates are independent
sigma = 1 #variance for error term of the outcome model

#empty data structures

Delta_b = NA
outcomes = matrix(NA,n,M)
treatments = matrix(NA,n,M)

############################### continuousdata ##########################################
continuousdata = list()
set.seed(123)

##iterate to create each dataset
for(m in 1:M){

  ###draw covariates from a multivariate normal distribution
  X = matrix( mvrnorm(n,mu=rep(0,p),Sigma=D), ncol = p)
  X = cbind(rep(1,n),X) #intercept
  
  #Treatment drawn based on a Bernoulli with probability = 0.5
  treatmentmatrix = matrix(rbinom(n,1,0.5),nrow=n,ncol=1)
  treatments[,m] = treatmentmatrix
  
  ###continuous outcome simulation
  error = rnorm(n,0,sigma) #random error term for outcome
  outcome = X%*%t(beta0) + X%*%t(beta1)*treatmentmatrix + error
  outcomes[,m] = outcome
  
  #ATE in benefit group
  G = data.frame(outcome=outcome,treat=treatmentmatrix,benefit=X%*%t(beta1))
  BG = G[G$benefit>0,]
  Delta_b = mean(BG[BG$treat==1,]$outcome)-mean(BG[BG$treat==0,]$outcome)
  #ATE in non benefit group
  #Delta_nb = outcome*(1-treatmentmatrix)
  
  #save simulated data
  continuousdata[[m]] = list()
  X = X[,-1]
  continuousdata[[m]][["covariates"]] = X
  continuousdata[[m]][["outcome"]] = outcome
  continuousdata[[m]][["treatment"]] = treatmentmatrix
  continuousdata[[m]][["benefitATE"]] = Delta_b
}

