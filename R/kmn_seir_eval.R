setwd("C:/Marci/CEU/ThesisDONTUSETHIS/CODE/R")
source("kmn_seir_functions.R")
library("optimParallel")

kmn_seir_error <- function(parameter_values, cases=flu$cases) {
  parameter_values = c(parameter_values, N = 504)
  pred <- kmn_seir(parameter_values, initial_values, 1:length(cases))
  return( sum( (pred$I - cases) ^ 2))
}
a = kmn_seir_error(parameter_values)
write.csv(a, "a.csv")
cases = flu$cases
kmn_optim <- function(parameter_values) {
  kmn_seir_equations <- function(time, variables, parameters) {
    
    with(as.list(c(variables, parameters)), {
      if (time <= tau_q+1){
        lagged = rep(0, 2)} #this is for S(t-tau_q) I_S(t-tau_q)
      else {
        lagged = deSolve :: lagvalue(time - tau_q)[c(1, 11)]
      }
      dS <- -(k*b*I+q*k*(1-b)*I_S)*S/N + (q*k*(1-b) * lagged[1] * lagged[2] )/N 
      
      dS_Q <- (q*k*(1-b)*S*I_S)/N - (q*k*(1-b) * lagged[1] * lagged[2] )/N
      
      dE1 <- (k*b*(I - q*I_S))*S/N - m*sigma*E1
      
      dE2 <- m * sigma * E1 - m * sigma * E2
      
      dI_A1 = m*sigma*E2 - n*gamma*I_A1 - P_I1
      
      dI_A2 <- n*gamma* I_A1 - n*gamma*I_A2 - P_I2
      
      dI_A3 <- n*gamma* I_A2 - n*gamma*I_A3 - P_I3
      
      dI_A <- dI_A1 + dI_A2 + dI_A3
      
      dI_S1 <- P_I1 - (n*gamma + d_I) * I_S1
      
      dI_S2 <- P_I2 + n * gamma * I_S1 - (n*gamma+d_I) * I_S2
      
      dI_S3 <- P_I3 + n * gamma * I_S3 - (n*gamma+d_I) * I_S3
      
      dI_S <- dI_S1+dI_S2+dI_S3
      
      dI = dI_A + dI_S
      
      dQ <- (q*k*b*S*I_S)/N + d_I*I_S
      
      dR <- n*gamma*(I_A3 + I_S3)
      
      dP_I1 = 0
      
      dP_I2 = 0
      
      dP_I3 = 0
      
      return(list(c(dS, dS_Q, dE1, dE2, dI_A1,dI_A2,dI_A3, dI_S1,dI_S2,dI_S3, dI_S, dI_A, dI, dP_I1, dP_I2, dP_I3, dQ, dR)))
    })
  }
  initial_values <- c(
    S=762,
    S_Q=0,
    E1 = 0,
    E2 = 0,
    I_A1 = 0,
    I_A2 = 0,
    I_A3 = 0,
    I_S1 = 1,
    I_S2 = 0,
    I_S3 = 0,
    I_S = 1,
    I_A =0,
    I=1,
    P_I1 = 0.1,
    P_I2 = 0.2,
    P_I3 = 0.14,
    Q=0,
    R=0
  )
  cases=read.csv("flu.csv")$cases
  out = deSolve :: dede(
    y=initial_values,
    times=1:length(cases),
    func=kmn_seir_equations,
    parms = c(parameter_values, N=763)
  )
  out = as.data.frame(out)
  return(sum((out$I - cases)^2))
}
kmn_optim(parameter_values)
###for flu

parameter_values <- c(
  k=100,
  b=0.04,
  q=0.02,
  sigma=1,
  m=2,
  gamma=1/5,
  n=3,
  d_I=0.5,
  tau_q = 5
)



initial_params = parameter_values
cl <- makeCluster(detectCores())
start = proc.time()
# the following runs for approx 1.5 min for the flu dataset
optimum <- optimParallel(par=initial_params, fn=kmn_optim, parallel=list(cl=cl))
end = proc.time()

end - start

best_params <- optimum$par
pred = kmn_seir(c(best_params, N=763), initial_values, 1:length(cases))
with(flu, plot(day, cases, pch = 19, col = "red", ylim = c(0, 600)))
# the model-predicted prevalences:
with(pred, lines(time, I, col = "red", type = "o"))
# the "errors":
segments(1:length(cases), flu$cases, pred$time, pred$I)
