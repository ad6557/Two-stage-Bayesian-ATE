#load packages + required functions
require(MASS)
seeed=123
set.seed(seeed)
###########
p = 10
beta0 = matrix(c(1,5,-15,15,-5,10,1,2,3,-4,-1),nrow=1)
D = diag(rep(1,p)) #correlation matrix for covaritates: all covariates are independent

############################## generate true ATE ###############################
N=10000000
###draw covariates from a multivariate normal distribution
X = matrix( mvrnorm(N,mu=rep(0,p),Sigma=D), ncol = p)
X = cbind(rep(1,N),X) #intercept

#Treatment drawn based on a Bernoulli with probability = 0.5
treatmentmatrix = matrix(c(rep(0,N/2),rep(1,N/2)),nrow=N,ncol=1)

Y = X%*%t(beta0) + 
  (35*sin(pi*X[,1]*X[,2])-3*(X[,3]-0.5)^2+3*X[,4]^2-X[,5] + 
     (sin(pi*X[,6]*X[,7])-(X[,8]+0.1)^2+X[,9]^2-0.5*X[,10]))*treatmentmatrix
Y1 = X%*%t(beta0) + 
  (35*sin(pi*X[,1]*X[,2])-3*(X[,3]-0.5)^2+3*X[,4]^2-X[,5] + 
     (sin(pi*X[,6]*X[,7])-(X[,8]+0.1)^2+X[,9]^2-0.5*X[,10]))*1 # if all treated
# hist(Y1)

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

n = 400 #sample size
M = 200 #replications

for(m in 1:M){
  
  ###draw covariates from a multivariate normal distribution
  X = matrix( mvrnorm(n,mu=rep(0,p),Sigma=D), ncol = p)
  X = cbind(rep(1,n),X) #intercept
  
  #Treatment drawn based on a Bernoulli with probability = 0.5
  treatment =matrix(c(rep(0,n/2),rep(1,n/2)),nrow=n,ncol=1)
  
  ### complex
  error = rnorm(n,0,1)
  Y = X%*%t(beta0) + 
    (35*sin(pi*X[,1]*X[,2])-3*(X[,3]-0.5)^2+3*X[,4]^2-X[,5] + 
       (sin(pi*X[,6]*X[,7])-(X[,8]+0.1)^2+X[,9]^2-0.5*X[,10]))*treatment + error
  
  outcome  <- Y
  
  #save simulated data
  continuousdata[[m]] = list()
  X = as.matrix(X[,-1],n,p)
  continuousdata[[m]][["covariates"]] = X[,1:5]
  continuousdata[[m]][["outcome"]] = outcome
  continuousdata[[m]][["treatment"]] = treatment
}


