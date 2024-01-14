#load packages + required functions
require(MASS)
seeed=123
set.seed(seeed)
###########
p = 5
beta0 = matrix(c(0,1,0.5,0,-0.5,-1),nrow=1)
beta1 = matrix(c(0,-0.5,0.5,-1,1,0),nrow=1)
D = diag(rep(1,p)) #correlation matrix for covaritates: all covariates are independent

############################## generate true ATE ###############################
N=100000000
###draw covariates from a multivariate normal distribution
X = matrix( mvrnorm(N,mu=rep(0,p),Sigma=D), ncol = p)
X = cbind(rep(1,N),X) #intercept

#Treatment drawn based on a Bernoulli with probability = 0.5
treatmentmatrix = matrix(c(rep(0,N/2),rep(1,N/2)),nrow=N,ncol=1)

###continuous outcome simulation
LP = X%*%t(beta0) + X%*%t(beta1)*treatmentmatrix
hist(pnorm(LP))
# logit(pi) = log(pi/(1-pi)) = x*b
# pi = 1/(1+exp(-x*b))
# plogis(x*b) = 1/(1+exp(-x*b))
outcome  <- rbinom(N, 1, pnorm(LP))


G = data.frame(outcome=outcome,treat=treatmentmatrix)
BG = G[LP>0,]
p1.benefit <- mean(BG[which(BG$treat==1),]$outcome)
p0.benefit <- mean(BG[which(BG$treat==0),]$outcome)

# benefit.log.OR <- log((p1.benefit/(1-p1.benefit))/(p0.benefit/(1-p0.benefit)))
PATE = p1.benefit - p0.benefit

############################### binarydata ##########################################
binarydata = list()

n = 200 #sample size
M = 500 #replications

for(m in 1:M){
  
  ###draw covariates from a multivariate normal distribution
  X = matrix( mvrnorm(n,mu=rep(0,p),Sigma=D), ncol = p)
  X = cbind(rep(1,n),X) #intercept
  
  #Treatment drawn based on a Bernoulli with probability = 0.5
  treatment = matrix(rbinom(n,1,0.5),nrow=n,ncol=1)
  
  ###continuous outcome simulation
  LP = X%*%t(beta0) + X%*%t(beta1)*treatment
  outcome  <- rbinom(n, 1, plogis(LP))
  
  #save simulated data
  binarydata[[m]] = list()
  X = as.matrix(X[,-1],n,p)
  binarydata[[m]][["covariates"]] = X
  binarydata[[m]][["outcome"]] = outcome
  binarydata[[m]][["treatment"]] = treatment
  binarydata[[m]][["true"]] = PATE
}


