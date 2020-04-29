setwd("C:/Marci/CEU/ThesisDONTUSETHIS/CODE/R")
source("MN-SEIR-MODELS/m3n5.R")
flu = read.csv("flu.csv")
cases = flu$cases

parameter_values <- c(
  k=100,
  b=0.05,
  q=0.02,
  sigma=1,
  gamma=1/5,
  d_I=0.5,
  tau_d = 2
)


initial_params = parameter_values
cl <- makeCluster(detectCores())
start = proc.time()
optimum <- optimParallel(par=initial_params,
                         fn=mn_optim,
                         lower=rep(0, 7),
                         control=list(maxit=150, trace=6, fnscale=1),
                         parallel=list(cl=cl))
end = proc.time()
stopCluster(cl)

end - start

