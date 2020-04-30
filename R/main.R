setwd("C:/Marci/CEU/ThesisDONTUSETHIS/CODE/R")
library(MCS)
library(optimParallel)
path <- function(m,n) {
  return(paste0("MN-SEIR-MODELS/m", m, "n", n, ".R"))
}
parscale.parameters <- function(par, scale, fix = 1){
  #check if length of scale is equal to par
  if(length(par) !=  length(scale)){
    stop("parscale.parameters has parameter and scaling vectors of different sizes.")
  }
  
  if(any(scale == 0)){
    scale[scale == 0] <- 1
  }
  #parscale fixes the larged par/parscale value to deviate only 10 percent, others can then vary
  #get fixed and maximal value
  fix.value <- par[fix]
  if(fix.value == 0){
    fix.value <- 1
  }
  
  max.value <- abs(par[which.max(abs(par))[1]])
  if(max.value == 0){
    print("Warning: Vector contains only zeroes. Scaling is set to a default of 1.")
    max.value <- 1
  }
  
  #fill in scaling vector
  par.scale <- scale/(0.1 * fix.value * max.value)
  par.scale[fix] <- 1 / max.value
  
  return(abs(par.scale))
}

parameter_values <- c(
  min_contract_size=10,
  lambda=10,
  b=0.08,
  q=0.9,
  sigma=1/5,
  gamma=1/14,
  tau_d=1,
  d_I=1
)
scale=c(1e-3, 0.05, 0.05, 0.1, 0.1, 0.2)
p.scale=parscale.parameters(parameter_values, scale)
vec=c(min_contract)

pars <- vector(mode="list", length=25)
vec = c()
for (m in 1:5) {
  for (n in 1:5) {
    vec <- c(vec, paste0("m",m,"n",n))
  }
}
names(pars) <- vec


for (m in 1:5) {
  for (n in 1:5) {
    source(path(m,n))
    initial_params = parameter_values
    cl <- makeCluster(detectCores())
    optimum <- optimParallel(par=initial_params,
                             fn=mn_optim,
                             lower=c(0.1,0.1, 0.001, 0.001, 0.05, 0.05, 0.001, 0.001),
                             control=list(maxit=150, trace=6, fnscale=25000000),
                             parallel=list(cl=cl))
    residuals=mn_pred(optimum$par)
    
    pars[[paste0("m",m,"n",n)]] = residuals
    print(paste("residuals for m=",m,", n=",n, "is ", residuals))
    stopCluster(cl)
    
  }
}
Loss = as.data.frame(pars)
MCS <- MCSprocedure(Loss=Loss,alpha=0.2,B=5000,statistic='Tmax',cl=cl)

