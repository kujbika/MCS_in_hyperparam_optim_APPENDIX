setwd("C:/Marci/CEU/ThesisDONTUSETHIS/CODE/R")
library(MCS)
path <- function(m,n) {
  return(paste0("MN-SEIR-MODELS/m", m, "n", n, ".R"))
}


parameter_values <- c(
  k=350,
  b=0.05,
  q=0.02,
  sigma=1,
  gamma=1/5,
  d_I=0.5,
  tau_d = 2
)


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
                             lower=rep(0, 7),
                             upper = c(700, 1, 1, Inf, Inf, 1, Inf),
                             control=list(maxit=150, trace=T, fnscale=1),
                             parallel=list(cl=cl))
    residuals=mn_pred(optimum$par)
    pars[[paste0("m",m,"n",n)]] = residuals
    print(paste("residuals for m=",m,", n=",n, "is ", residuals))
    stopCluster(cl)
    
  }
}
Loss = as.data.frame(pars)
MCS <- MCSprocedure(Loss=Loss,alpha=0.2,B=5000,statistic='Tmax',cl=cl)

