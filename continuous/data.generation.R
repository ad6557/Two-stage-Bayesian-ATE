require(MASS)

p = 10
beta0 = matrix(c(0, 0.5, 1, 0.75, 1, 0.5,-0.5,-1,-0.75,-1,-0.5),nrow=1)
beta1 = matrix(c(0, 0.5,-0.5,0.65,-0.65,1,-1,0.75,-0.75,0.25,-0.25),nrow=1)
D = diag(rep(1,p)) #correlation matrix for covaritates

############################## generate true ATE ###############################
N=10000000
###draw covariates from a multivariate normal distribution
X = matrix( mvrnorm(N,mu=rep(0,p),Sigma=D), ncol = p)
X = cbind(rep(1,N),X)

#Treatment drawn based on a Bernoulli with probability = 0.5
treatmentmatrix = matrix(c(rep(0,N/2),rep(1,N/2)),nrow=N,ncol=1)

### linear
Y = X%*%t(beta0) + X%*%t(beta1)*treatmentmatrix
Y1 = X%*%t(beta0) + X%*%t(beta1)*1 # if all treated

minCond = 0
G = data.frame(Y=Y,TRT=treatmentmatrix)
NLR = G[Y1<=minCond,]
LR = G[Y1>minCond,]

# didn't use "bzGetSubgrpRaw" because it is too slow when data size large
ATE_all = mean(G[G$TRT==1,]$Y)-mean(G[G$TRT==0,]$Y)
ATE_NLR = mean(NLR[NLR$TRT==1,]$Y)-mean(NLR[NLR$TRT==0,]$Y)
ATE_LR = mean(LR[LR$TRT==1,]$Y)-mean(LR[LR$TRT==0,]$Y)

############################### continuousdata ##########################################
continuousdata = list()

n = 500 #sample size
M = 200 #replications

for(m in 1:M){
  
  # covariate generation
  X = matrix( mvrnorm(n,mu=rep(0,p),Sigma=D), ncol = p)
  X = cbind(rep(1,n),X) #intercept
  # treatment generation
  treatment =matrix(c(rep(0,n/2),rep(1,n/2)),nrow=n,ncol=1)
  # outcome generation
  error = rnorm(n,0,1)
  outcome = X%*%t(beta0) + X%*%t(beta1)*treatment + error
  
  #save simulated data
  continuousdata[[m]] = list()
  covariates = as.matrix(X[,-1],n,p)
  continuousdata[[m]][["covariates"]] = covariates
  continuousdata[[m]][["outcome"]] = outcome
  continuousdata[[m]][["treatment"]] = treatment
}


