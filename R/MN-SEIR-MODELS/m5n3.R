library("optimParallel")

mn_optim <- function(parameter_values) {
  factorial <- function(n, acc=1) {
    if (n <= 0) {
      return(acc)
    } else {
      return(factorial(n-1, acc * n))
    }
  }
  P.gen <- function(m, sigma, E_m_lagged, tau_d, n, gamma, i) {
    return(m * sigma * E_m_lagged * exp(-n*gamma*tau_d)*((n*gamma*tau_d)^(i-1))/factorial(i-1))
  }
  
  mn_seir_equations <- function(time, variables, parameters) {
    
    with(as.list(c(variables, parameters)), {
      if (time <= tau_q+1){
        lagged = rep(0, 2)} #this is for S(t-tau_q) and I_S(t-tauq)
      else {
        lagged = deSolve :: lagvalue(time - tau_q)[c(1, 14)]
      }
      dS <- -(k*b*I+q*k*(1-b)*I_S)*S/N + (q*k*(1-b) * lagged[1] * lagged[2] )/N 
      
      dS_Q <- (q*k*(1-b)*S*I_S)/N - (q*k*(1-b) * lagged[1] * lagged[2] )/N
      
      dE1 <- (k*b*(I - q*I_S))*S/N - m*sigma*E1
      
      dE2 <- m * sigma * (E1 - E2)
      
      dE3 <- m * sigma * (E2 - E3)
      
      dE4 <- m * sigma * (E3 - E4)
      
      dE5 <- m * sigma * (E4 - E5)
      
      
      if (time <= tau_d+1){
        lagged1 = rep(0,2)} #this is for E_m(t-tau_d)
      else {
        lagged1 = deSolve :: lagvalue(time - tau_d)[c(6,7)]
      }
      
      dP_I1 = P.gen(m, sigma, m * sigma * (lagged1[1] - lagged1[2]), tau_d, n, gamma, 1)
      
      dP_I2 = P.gen(m, sigma, m * sigma * (lagged1[1] - lagged1[2]), tau_d, n, gamma, 2)
      
      dP_I3 = P.gen(m, sigma, m * sigma * (lagged1[1] - lagged1[2]), tau_d, n, gamma, 3)
      
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
      
      return(list(c(dS, dS_Q, dE1, dE2, dE3, dE4, dE5, dI_A1,dI_A2,dI_A3, dI_S1,dI_S2,dI_S3, dI_S, dI_A, dI, dP_I1, dP_I2, dP_I3, dQ, dR)))
    })
  }
  initial_values <- c(
    S=762,
    S_Q=0,
    E1 = 0,
    E2 = 0,
    E3=0,
    E4=0,
    E5=0,
    I_A1 = 0,
    I_A2 = 0,
    I_A3 = 0,
    I_S1 = 1,
    I_S2 = 0,
    I_S3 = 0,
    I_S = 1,
    I_A =0,
    I=1,
    P_I1 = 0.0,
    P_I2 = 0.0,
    P_I3 = 0.0,
    Q=0,
    R=0
  )
  cases=read.csv("flu.csv")$cases
  out = deSolve :: dede(
    y=initial_values,
    times=1:length(cases),
    func=mn_seir_equations,
    parms = c(parameter_values, N=763, m=5, n=3, tau_q=5),
    method = "impAdams",
    control = list(interpol=2)
  )
  out = as.data.frame(out)
  return(sum(((out$I + out$Q) - cases)^2))
}

