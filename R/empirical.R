library("deSolve")
library(lazyeval)
caller <- function(str, idx) {
  return (get(paste0(str,idx)))
}

LazyEval_cheatcode_dfe <- function(m,n) {
  str = paste()
}

kmn_seir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    if (time <= tau_q+1){
      lagged = rep(0, 2)} #this is for S(t-tau_q) I_S(t-tau_q)
    else {
      lagged = lagvalue(time - tau_q)[c(1, 11)]
    }
    dS <- -(k*b*I+q*k*(1-b)*I_S)*S/N + (q*k*(1-b) * lagged[1] * lagged[2] )/N 
    dS_Q <- (q*k*(1-b)*S*I_S)/N - (q*k*(1-b) * lagged[1] * lagged[2] )/N
    dE1 <- (k*b*(I - q*I_S))*S/N - m*sigma*E1
    dE2 <- m * sigma * E1 - m * sigma * E2
    #for (i in 2:m) {
    #  name <- paste0("dE", i)
    #  assign(name, m * sigma * caller("E", (i-1) ) - m * sigma * caller("E",i))
    #}
    #dE_all = unlist(lapply(1:m, function(x) caller("dE",x)))
    
    dI_A1 = m*sigma*E2 - n*gamma*I_A1 - P_I1
    
    dI_A2 <- n*gamma* I_A1 - n*gamma*I_A2 - P_I2
    
    dI_A3 <- n*gamma* I_A2 - n*gamma*I_A3 - P_I3
    
    dI_A <- dI_A1 + dI_A2 + dI_A3
    #for (i in 2:n) {
    #  name <- paste0("dI_A", i)
    #  assign(name, n*gamma* caller("I_A",i-1) - n*gamma* caller("I_A",i) - caller("P_I",i) )
    #}
    #dI_A_all = unlist(lapply(1:n, function(x) caller("dI_A",x)))
    
    #dI_A = sum(dI_A_all)
    
    dI_S1 <- P_I1 - (n*gamma + d_I) * I_S1
    dI_S2 <- P_I2 + n * gamma * I_S1 - (n*gamma+d_I) * I_S2
    dI_S3 <- P_I3 + n * gamma * I_S3 - (n*gamma+d_I) * I_S3
    dI_S <- dI_S1+dI_S2+dI_S3
    #for (i in 2:n) {
    #  name <- paste0("dI_S",i)
    #  assign(name, caller("P_I",i) + n * gamma * caller("I_S",i-1) - (n*gamma+d_I) * caller("I_S",i)) 
    #}
    #dI_S_all = unlist(lapply(1:n, function(x) caller("dI_S",x)))
    
    #dI_S = sum(dI_S_all)
    
    dI = dI_A + dI_S
    
    dQ <- (q*k*b*S*I_S)/N + d_I*I_S
    
    dR <- n*gamma*(I_A3 + I_S3)
    

    dP_I1 = 0
    dP_I2 = 0
    dP_I3 = 0

    
    return(list(c(dS, dS_Q, dE1, dE2, dI_A1,dI_A2,dI_A3, dI_S1,dI_S2,dI_S3, dI_S, dI_A, dI, dP_I1, dP_I2, dP_I3, dQ, dR)))
  })
}

parameter_values <- c(
  k=100,
  b=0.004,
  q=1,
  N=504,
  m=2,
  sigma=0.1,
  n=3,
  gamma=1/14,
  d_I=0.8,
  tau_q = 14
)

vector_maker <- function(str, MorN) {
  alls = unlist(lapply(1:parameter_values[MorN], function (x) assign(paste0(str,x), 0)))
  names_ = unlist(lapply(1:parameter_values[MorN], function (x) paste0(str,x)))
  names(alls) = names_
  return(alls)
}

p_gen <- function(m, sigma, E_m_lagged, tau_d, n, gamma, i) {
  return(m * sigma * E_m_lagged * exp(-n*gamma*tau_d)*((n*gamma*tau_d)^(i-1))/(i-1))
}

initial_values <- c(
  S=500,
  S_Q=0,
  E1 = 1,
  E2 = 1,
  I_A1 = 1,
  I_A2 = 0,
  I_A3 = 0,
  I_S1 = 1,
  I_S2 = 0,
  I_S3 = 0,
  I_S = 1,
  I_A =1,
  I=2,
  P_I1 = 0.0,
  P_I2 = 0.0,
  P_I3 = 0.0,
  Q=0,
  R=0
)

time_values <- seq(0:100)
res = dede(
  y=initial_values,
  times=time_values,
  func=kmn_seir_equations,
  parms = parameter_values
)

res = as.data.frame(res)
with(res, {
  # plotting the time series of susceptibles:
  plot(time, S, type = "l", col = "blue",
       xlab = "time (days)", ylab = "number of people")
  # adding the time series of infectious:
  lines(time, I, col = "red")
  # adding the time series of recovered:
  lines(time, R, col = "green")
  lines(time, Q, col="black")
  lines(time, S_Q, col="purple")
})

# adding a legend:
legend("right", c("susceptibles", "infectious", "recovered", "quarantine"),
       col = c("blue", "red", "green", "black"), lty = 1, bty = "n")

write.csv(res, "C:/OwnWork/thesis/results.csv")


kmn_seir <- function(parameter_values, initial_values, time_values) {
  res = dede(
    y=initial_values,
    times=time_values,
    func=kmn_seir_equations,
    parms = parameter_values
  )
  return(as.data.frame(res))
}

flu <- read.table("https://bit.ly/2vDqAYN", header = TRUE)
flu
preds <- kmn_seir(parameter_values, initial_values, flu$day)









