setwd("C:/Marci/CEU/ThesisDONTUSETHIS/CODE/R")

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
res <- data.frame(matrix(rep(0, 25), nrow=5))
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
    res[m, n] = optimum$value
    print(paste0("m is ",m,", n is ",n))
    print(paste("rmse is ", optimum$value))
    stopCluster(cl)
    
  }
}

