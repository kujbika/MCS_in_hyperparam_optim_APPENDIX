setwd("C:/Marci/CEU/ThesisDONTUSETHIS/CODE/R/Covid19_Prediction")
load("covidPred_secondoptim.RData")
library(MCS)
library(optimParallel)
path <- function(m,n) {
  return(paste0("MN-SEIR-MODELS/m", m, "n", n, ".R"))
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


pars <- vector(mode="list", length=25)
vec = c()
for (m in 1:5) {
  for (n in 1:5) {
    vec <- c(vec, paste0("m",m,"n",n))
  }
}
names(pars) <- vec

coeffs <- vector(mode="list", length=25)
names(coeffs) <- vec



for (m in 1:5) {
  for (n in 1:5) {
    source(path(m,n))
    initial_params = parameter_values
    cl <- makeCluster(detectCores())
    optimum <- optimParallel(par=initial_params,
                             fn=mn_optim,
                             lower=rep(0,7),
                             control=list(maxit=150, trace=6, fnscale=250000000),
                             parallel=list(cl=cl))
    residuals=mn_pred(optimum$par)
    
    pars[[paste0("m",m,"n",n)]] = residuals
    coeffs[[paste0("m",m,"n",n)]] = optimum$par
    print(paste("residuals for m=",m,", n=",n, "is ", residuals))
    stopCluster(cl)
    
  }
}
Loss = as.data.frame(pars)
MCS <- MCSprocedure(Loss=Loss,alpha=0.2,B=5000,statistic='Tmax',cl=cl)

