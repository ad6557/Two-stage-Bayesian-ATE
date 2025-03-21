require(MASS)

p = 10
beta0 = matrix(c(0, 0.5, 1, 0.75, 1, 0.5,-0.5,-1,-0.75,-1,-0.5),nrow=1)
beta1 = matrix(c(0, 0.5,-0.5,0.65,-0.65,0.35,-0.35,0.75,-0.75,0.25,-0.25),nrow=1)
D = diag(rep(1,p)) #correlation matrix for covaritates

############################## generate true ATE ###############################
N=10000000
###draw covariates from a multivariate normal distribution
X = matrix( mvrnorm(N,mu=rep(0,p),Sigma=D), ncol = p)
X = cbind(rep(1,N),X)

#Treatment drawn based on a Bernoulli with probability = 0.5
treatmentmatrix = matrix(c(rep(0,N/2),rep(1,N/2)),nrow=N,ncol=1)

### linear
LP = X%*%t(beta0) + X%*%t(beta1)*treatmentmatrix
LP1 = X%*%t(beta0) + X%*%t(beta1)*1 # if all treated
outcome  <- rbinom(N, 1, plogis(LP))

minCond = 0
G = data.frame(Y=outcome,TRT=treatmentmatrix)
NLR = G[LP1<=minCond,]
LR = G[LP1>minCond,]

OR_all = fisher.test(table(G))[["estimate"]][["odds ratio"]]
OR_NLR = fisher.test(table(NLR))[["estimate"]][["odds ratio"]]
OR_LR = fisher.test(table(LR))[["estimate"]][["odds ratio"]]
log_NLR=log(OR_NLR);log_LR=log(OR_LR);log_all=log(OR_all);

############################### binarydata ##########################################
binarydata = list()

n = 500 #sample size
M = 200 #replications

for(m in 1:M){
  
  # covariate generation
  X = matrix( mvrnorm(n,mu=rep(0,p),Sigma=D), ncol = p)
  X = cbind(rep(1,n),X) #intercept
  # treatment generation
  treatment =matrix(c(rep(0,n/2),rep(1,n/2)),nrow=n,ncol=1)
  # binary outcome generation
  LP = X%*%t(beta0) + X%*%t(beta1)*treatment
  outcome  <- rbinom(n, 1, plogis(LP))
  
  #save simulated data
  binarydata[[m]] = list()
  covariates = as.matrix(X[,-1],n,p)
  binarydata[[m]][["covariates"]] = covariates
  binarydata[[m]][["outcome"]] = outcome
  binarydata[[m]][["treatment"]] = treatment
}


