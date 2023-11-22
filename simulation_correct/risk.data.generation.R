#load packages + required functions
require(MASS)

###########
p = 5
beta0 = matrix(c(0.1,seq(0.2,1,length=p)),nrow=1) #covariates for PS model
#beta1 = matrix(c(0.3,c(-1.1,-0.8,-0.5,-0.2,0.1,0.4,0.7,1,1.3)),nrow=1)
beta1 = matrix(c(0.3,c(-0.8,-0.5,-0.1,0.5,1)),nrow=1)
#beta1 = matrix(c(0.8,c(-0.3,0.3)),nrow=1)
D = diag(rep(1,p)) #correlation matrix for covaritates: all covariates are independent
sigma = 1 #variance for error term of the outcome model


############################## generate true ATE ###############################
N=100000000
###draw covariates from a multivariate normal distribution
X = matrix( mvrnorm(N,mu=rep(0,p),Sigma=D), ncol = p)
X = cbind(rep(1,N),X) #intercept

#Treatment drawn based on a Bernoulli with probability = 0.5
treatmentmatrix = matrix(rbinom(N,1,0.5),nrow=N,ncol=1)

###continuous outcome simulation
outcome = X%*%t(beta0) + X%*%t(beta1)*treatmentmatrix #+ error


G = data.frame(outcome=outcome,treat=treatmentmatrix,benefit=X%*%t(beta1))
#ATE in benefit group
BG = G[G$benefit>0,]
Delta_b = mean(BG[BG$treat==1,]$outcome)-mean(BG[BG$treat==0,]$outcome)
#ATE in non benefit group
#NBG = G[G$benefit<=0,]
#Delta_nb = mean(NBG[NBG$treat==1,]$outcome)-mean(NBG[NBG$treat==0,]$outcome)


############################### continuousdata ##########################################
continuousdata = list()
set.seed(1234)

n =100 #sample size
M = 200 #replications
outcomes = matrix(NA,n,M)
treatments = matrix(NA,n,M)

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
  
  
  #save simulated data
  continuousdata[[m]] = list()
  X = X[,-1]
  continuousdata[[m]][["covariates"]] = X
  continuousdata[[m]][["outcome"]] = outcome
  continuousdata[[m]][["treatment"]] = treatmentmatrix
  continuousdata[[m]][["benefitATE"]] = Delta_b
}



