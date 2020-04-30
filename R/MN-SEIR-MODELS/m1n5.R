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
        lagged = deSolve :: lagvalue(time - tau_q)[c(1, 2)]
      }
      dS <- -(K*b*I+q*K*(1-b)*I_S)*S/N + (q*K*(1-b) * lagged[1] * lagged[2] )/N 
      
      dS_Q <- (q*K*(1-b)*S*I_S)/N - (q*K*(1-b) * lagged[1] * lagged[2] )/N
      
      dE1 <- (K*b*(I - q*I_S))*S/N - m*sigma*E1
      
      if (time <= tau_d+1){
        lagged1 = rep(0,4)
        dEm = 0} #this is for E_m(t-tau_d)
      else {
        lagged1 = deSolve :: lagvalue(time - tau_d)[c(1,2,3,5)]
        dEm <- (1/N)*K*b*(lagged1[3] - q*lagged1[2]) * lagged1[1] - m*sigma*lagged1[4]
      }
      
      dP_I1 = P.gen(m, sigma, dEm, tau_d, n, gamma, 1)
      
      dP_I2 = P.gen(m, sigma, dEm, tau_d, n, gamma, 2)
      
      dP_I3 = P.gen(m, sigma, dEm, tau_d, n, gamma, 3)
      
      dP_I4 = P.gen(m, sigma, dEm, tau_d, n, gamma, 4)
      
      dP_I5 = P.gen(m, sigma, dEm, tau_d, n, gamma, 5)
      
      
      dI_A1 = m*sigma*E1 - n*gamma*I_A1 - P_I1
      
      dI_A2 <- n*gamma* I_A1 - n*gamma*I_A2 - P_I2 
      
      dI_A3 <- n*gamma* I_A2 - n*gamma*I_A3 - P_I3
      
      dI_A4 <- n*gamma* I_A3 - n*gamma*I_A4 - P_I4
      
      dI_A5 <- n*gamma* I_A4 - n*gamma*I_A5 - P_I5
      
      
      dI_A <- dI_A1 + dI_A2 + dI_A3 + dI_A4 + dI_A5
      
      dI_S1 <- P_I1 - (n*gamma + d_I) * I_S1
      
      dI_S2 <- P_I2 + n * gamma * I_S1 - (n*gamma+d_I) * I_S2
      
      dI_S3 <- P_I3 + n * gamma * I_S2 - (n*gamma+d_I) * I_S3
      
      dI_S4 <- P_I4 + n * gamma * I_S3 - (n*gamma+d_I) * I_S4
      
      dI_S5 <- P_I5 + n * gamma * I_S4 - (n*gamma+d_I) * I_S5
      
      dI_S <- dI_S1 + dI_S2 + dI_S3 + dI_S4 + dI_S5 
      
      dI <- dI_A + dI_S
      
      dQ <- (q*K*b*S*I_S)/N + d_I*I_S
      
      dR_A <- n*gamma*I_A5
      dR <- n*gamma*I_S5
      dK <- -(K-min_contract_size)/lambda
      return(list(c(dS, dI_S, dI, dS_Q, dE1, dI_A1, dI_A2, dI_A3, dI_A4, dI_A5, dI_S1, dI_S2, dI_S3, dI_S4, dI_S5, dI_A, dP_I1, dP_I2, dP_I3, dP_I4, dP_I5, dQ, dR, dK, dR_A)))
    })
  }
  initial_values <- c(
    S=19449999,
    I_S = 1,
    I=1,
    S_Q=0,
    E1 = 0,
    I_A1 = 0,
    I_A2 = 0,
    I_A3 = 0,
    I_A4 = 0,
    I_A5 = 0,
    I_S1 = 1,
    I_S2 = 0,
    I_S3 = 0,
    I_S4 = 0,
    I_S5 = 0,
    I_A =0,
    P_I1 = 0.0,
    P_I2 = 0.0,
    P_I3 = 0.0,
    P_I4 = 0.0,
    P_I5 = 0.0,
    Q=0,
    R=0,
    K=50,
    R_A=0
  )
  cases=read.csv("curve.csv")$x
  cases=head(cases, length(cases) * 0.75)
  out = deSolve :: dede(
    y=initial_values,
    times=1:length(cases),
    func=mn_seir_equations,
    parms = c(parameter_values, N=19450000, m=1, n=5, tau_q=14),
    method = "impAdams",
    control = list(interpol=2)
  )
  out = as.data.frame(out)
  return(sum(((out$I_S + out$Q + out$R) - cases)^2) / length(cases))
}

mn_pred <- function(parameter_values) {
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
        lagged = deSolve :: lagvalue(time - tau_q)[c(1, 2)]
      }
      dS <- -(K*b*I+q*K*(1-b)*I_S)*S/N + (q*K*(1-b) * lagged[1] * lagged[2] )/N 
      
      dS_Q <- (q*K*(1-b)*S*I_S)/N - (q*K*(1-b) * lagged[1] * lagged[2] )/N
      
      dE1 <- (K*b*(I - q*I_S))*S/N - m*sigma*E1
      
      if (time <= tau_d+1){
        lagged1 = rep(0,4)
        dEm = 0} #this is for E_m(t-tau_d)
      else {
        lagged1 = deSolve :: lagvalue(time - tau_d)[c(1,2,3,5)]
        dEm <- (1/N)*K*b*(lagged1[3] - q*lagged1[2]) * lagged1[1] - m*sigma*lagged1[4]
      }
      
      dP_I1 = P.gen(m, sigma, dEm, tau_d, n, gamma, 1)
      
      dP_I2 = P.gen(m, sigma, dEm, tau_d, n, gamma, 2)
      
      dP_I3 = P.gen(m, sigma, dEm, tau_d, n, gamma, 3)
      
      dP_I4 = P.gen(m, sigma, dEm, tau_d, n, gamma, 4)
      
      dP_I5 = P.gen(m, sigma, dEm, tau_d, n, gamma, 5)
      
      
      dI_A1 = m*sigma*E1 - n*gamma*I_A1 - P_I1
      
      dI_A2 <- n*gamma* I_A1 - n*gamma*I_A2 - P_I2 
      
      dI_A3 <- n*gamma* I_A2 - n*gamma*I_A3 - P_I3
      
      dI_A4 <- n*gamma* I_A3 - n*gamma*I_A4 - P_I4
      
      dI_A5 <- n*gamma* I_A4 - n*gamma*I_A5 - P_I5
      
      
      dI_A <- dI_A1 + dI_A2 + dI_A3 + dI_A4 + dI_A5
      
      dI_S1 <- P_I1 - (n*gamma + d_I) * I_S1
      
      dI_S2 <- P_I2 + n * gamma * I_S1 - (n*gamma+d_I) * I_S2
      
      dI_S3 <- P_I3 + n * gamma * I_S2 - (n*gamma+d_I) * I_S3
      
      dI_S4 <- P_I4 + n * gamma * I_S3 - (n*gamma+d_I) * I_S4
      
      dI_S5 <- P_I5 + n * gamma * I_S4 - (n*gamma+d_I) * I_S5
      
      dI_S <- dI_S1 + dI_S2 + dI_S3 + dI_S4 + dI_S5 
      
      dI <- dI_A + dI_S
      
      dQ <- (q*K*b*S*I_S)/N + d_I*I_S
      
      dR_A <- n*gamma*I_A5
      dR <- n*gamma*I_S5
      
      dK <- -(K-min_contract_size)/lambda
      return(list(c(dS, dI_S, dI, dS_Q, dE1, dI_A1, dI_A2, dI_A3, dI_A4, dI_A5, dI_S1, dI_S2, dI_S3, dI_S4, dI_S5, dI_A, dP_I1, dP_I2, dP_I3, dP_I4, dP_I5, dQ, dR, dK, dR_A)))
    })
  }
  initial_values <- c(
    S=19449999,
    I_S = 1,
    I=1,
    S_Q=0,
    E1 = 0,
    I_A1 = 0,
    I_A2 = 0,
    I_A3 = 0,
    I_A4 = 0,
    I_A5 = 0,
    I_S1 = 1,
    I_S2 = 0,
    I_S3 = 0,
    I_S4 = 0,
    I_S5 = 0,
    I_A =0,
    P_I1 = 0.0,
    P_I2 = 0.0,
    P_I3 = 0.0,
    P_I4 = 0.0,
    P_I5 = 0.0,
    Q=0,
    R=0,
    K=50,
    R_A=0
  )
  cases=read.csv("curve.csv")$x
  out = deSolve :: dede(
  y=initial_values,
  times=1:length(cases),
  func=mn_seir_equations,
  parms=c(parameter_values, N=19450000, m=1, n=5, tau_q=14),
  method="impAdams",
  control=list(interpol=2)
  )
  pred=tail(as.data.frame(out)[,c("I_S","Q", "R")], length(cases) * 0.25)
  residuals=apply(pred, 1, sum) - tail(cases, length(cases) * 0.25)
  return(residuals^2)
}
