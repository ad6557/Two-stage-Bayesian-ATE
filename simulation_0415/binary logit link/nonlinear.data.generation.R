#load packages + required functions
require(MASS)
seeed=123
set.seed(seeed)
###########
#p = 5
#beta0 = matrix(c(0,1,0.5,0,-0.5,-1),nrow=1)
p = 10
beta0 = matrix(c(0,1,0.5,0,-0.5,-1, 0.1,0.05,0,-0.05,-0.1),nrow=1)
D = diag(rep(0.16,p)) #correlation matrix for covaritates: all covariates are independent

############################## generate true ATE ###############################
N=10000000
###draw covariates from a multivariate normal distribution
X = matrix( mvrnorm(N,mu=rep(0,p),Sigma=D), ncol = p)
X = cbind(rep(1,N),X) #intercept

#Treatment drawn based on a Bernoulli with probability = 0.5
treatmentmatrix = matrix(c(rep(0,N/2),rep(1,N/2)),nrow=N,ncol=1)

### when p=5
# LP = X%*%t(beta0) + (sin(pi*X[,1]*X[,2])-(X[,3]+0.1)^2+X[,4]^2-X[,5]^3)*treatmentmatrix
# LP1 = X%*%t(beta0) + (sin(pi*X[,1]*X[,2])-(X[,3]+0.1)^2+X[,4]^2-X[,5]^3)*1 # if all treated

### when p=10
LP = X%*%t(beta0) + 
  (sin(pi*X[,1]*X[,2])-(X[,3]+0.1)^2+X[,4]^2-X[,5]^3 + 
     0.2*(sin(pi*X[,6]*X[,7])-(X[,8]+0.1)^2+X[,9]^2-X[,10]^3))*treatmentmatrix
LP1 = X%*%t(beta0) + 
  (sin(pi*X[,1]*X[,2])-(X[,3]+0.1)^2+X[,4]^2-X[,5]^3 + 
     0.2*(sin(pi*X[,6]*X[,7])-(X[,8]+0.1)^2+X[,9]^2-X[,10]^3))*1 # if all treated

hist(LP1)
# logit(pi) = log(pi/(1-pi)) = x*b
# pi = 1/(1+exp(-x*b))
# plogis(x*b) = 1/(1+exp(-x*b))
outcome  <- rbinom(N, 1, plogis(LP))

minCond = 0
G = data.frame(Y=outcome,TRT=treatmentmatrix)
NLR = G[LP1<=minCond,]
LR = G[LP1>minCond,]

# didn't use "bzGetSubgrpRaw" because it is too slow when data size large
OR_all = fisher.test(table(G))[["estimate"]][["odds ratio"]]
OR_NLR = fisher.test(table(NLR))[["estimate"]][["odds ratio"]]
OR_LR = fisher.test(table(LR))[["estimate"]][["odds ratio"]]
############################### binarydata ##########################################
binarydata = list()

n = 800 #sample size
M = 200 #replications

for(m in 1:M){
  
  ###draw covariates from a multivariate normal distribution
  X = matrix( mvrnorm(n,mu=rep(0,p),Sigma=D), ncol = p)
  X = cbind(rep(1,n),X) #intercept
  
  #Treatment drawn based on a Bernoulli with probability = 0.5
  treatment =matrix(c(rep(0,n/2),rep(1,n/2)),nrow=n,ncol=1)
  
  ### p=5
  #LP = X%*%t(beta0) + (sin(pi*X[,1]*X[,2])-(X[,3]+0.1)^2+X[,4]^2-X[,5]^3)*treatment
  ### p=10
  LP = X%*%t(beta0) + 
    (sin(pi*X[,1]*X[,2])-(X[,3]+0.1)^2+X[,4]^2-X[,5]^3 + 
       0.2*(sin(pi*X[,6]*X[,7])-(X[,8]+0.1)^2+X[,9]^2-X[,10]^3))*treatment
  
  outcome  <- rbinom(n, 1, plogis(LP))
  
  #save simulated data
  binarydata[[m]] = list()
  X = as.matrix(X[,-1],n,p)
  binarydata[[m]][["covariates"]] = X #[,1:5]
  binarydata[[m]][["outcome"]] = outcome
  binarydata[[m]][["treatment"]] = treatment
}


