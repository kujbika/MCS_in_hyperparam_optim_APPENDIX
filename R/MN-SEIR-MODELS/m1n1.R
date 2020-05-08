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
        lagged = rep(0, 3)} #this is for S(t-tau_q) and I_S(t-tauq)
      else {
        l = deSolve :: lagvalue(time - tau_q)
      lagged = c(l[c(1, 2)], rev(l)[3])
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
      dI_A1 = m*sigma*E1 - n*gamma*I_A1 - P_I1
      
      
      dI_A <- dI_A1
      
      dI_S1 <- P_I1 - (n*gamma + d_I) * I_S1
      
      dI_S <- dI_S1
      
      dI <- dI_A + dI_S
      
      dQ <- (q*K*b*S*I_S)/N + d_I*I_S
      
      dR <- n*gamma*(I_A1 + I_S1)
      dK <- -(K-min_contract_size)/lambda
      return(list(c(dS, dI_S, dI, dS_Q, dE1, dI_A1, dI_S1, dI_A, dP_I1, dQ, dR, dK)))
    })
  }
  initial_values <- c(
    S=9999998,
    I_S = 2,
    I=2,
    S_Q=0,
    E1 = 0,
    I_A1 = 0,
    I_S1 = 2,
    I_A =0,
    P_I1 = 0.0,
    Q=0,
    R=0,
    K=50
  )
  cases=read.csv("curve.csv")$x
  cases=head(cases, length(cases) * 0.75)
  out = deSolve :: dede(
    y=initial_values,
    times=1:length(cases),
    func=mn_seir_equations,
    parms = c(parameter_values, N=10000000, m=1, n=1, tau_q=14),
    method = "impAdams",
    control = list(interpol=2)
  )
  out = as.data.frame(out)
  return(sum(((out$I_S + out$Q + out$R - cases)^2) / length(cases)))
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
        lagged = rep(0, 3)} #this is for S(t-tau_q) and I_S(t-tauq)
      else {
        l = deSolve :: lagvalue(time - tau_q)
      lagged = c(l[c(1, 2)], rev(l)[3])
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
      dI_A1 = m*sigma*E1 - n*gamma*I_A1 - P_I1
      
      
      dI_A <- dI_A1
      
      dI_S1 <- P_I1 - (n*gamma + d_I) * I_S1
      
      dI_S <- dI_S1
      
      dI <- dI_A + dI_S
      
      dQ <- (q*K*b*S*I_S)/N + d_I*I_S
      
      dR <- n*gamma*(I_A1 + I_S1)
      dK <- -(K-min_contract_size)/lambda
      return(list(c(dS, dI_S, dI, dS_Q, dE1, dI_A1, dI_S1, dI_A, dP_I1, dQ, dR, dK)))
    })
  }
  initial_values <- c(
    S=9999998,
    I_S = 2,
    I=2,
    S_Q=0,
    E1 = 0,
    I_A1 = 0,
    I_S1 = 2,
    I_A =0,
    P_I1 = 0.0,
    Q=0,
    R=0,
    K=50
  )
  cases=read.csv("curve.csv")$x
  out = deSolve :: dede(
    y=initial_values,
    times=1:length(cases),
    func=mn_seir_equations,
    parms = c(parameter_values, N=10000000, m=1, n=1, tau_q=14),
    method = "impAdams",
    control = list(interpol=2)
  )
  pred = tail(as.data.frame(out)[,c("I_S", "Q", "R")], length(cases) * 0.25)
  residuals=apply(pred, 1, sum) - tail(cases, length(cases) * 0.25)
  return(residuals^2)
}

mn_dataframe <- function(parameter_values) {
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
        lagged = rep(0, 3)} #this is for S(t-tau_q) and I_S(t-tauq)
      else {
        l = deSolve :: lagvalue(time - tau_q)
        lagged = c(l[c(1, 2)], rev(l)[3])
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
      dI_A1 = m*sigma*E1 - n*gamma*I_A1 - P_I1
      
      
      dI_A <- dI_A1
      
      dI_S1 <- P_I1 - (n*gamma + d_I) * I_S1
      
      dI_S <- dI_S1
      
      dI <- dI_A + dI_S
      
      dQ <- (q*K*b*S*I_S)/N + d_I*I_S
      
      dR <- n*gamma*(I_A1 + I_S1)
      dK <- -(K-min_contract_size)/lambda
      return(list(c(dS, dI_S, dI, dS_Q, dE1, dI_A1, dI_S1, dI_A, dP_I1, dQ, dR, dK)))
    })
  }
  initial_values <- c(
    S=9999998,
    I_S = 2,
    I=2,
    S_Q=0,
    E1 = 0,
    I_A1 = 0,
    I_S1 = 2,
    I_A =0,
    P_I1 = 0.0,
    Q=0,
    R=0,
    K=50
  )
  cases=read.csv("curve.csv")$x
  cases=head(cases, length(cases) * 1)
  out = deSolve :: dede(
    y=initial_values,
    times=1:length(cases),
    func=mn_seir_equations,
    parms = c(parameter_values, N=10000000, m=1, n=1, tau_q=14),
    method = "impAdams",
    control = list(interpol=2)
  )
  out = as.data.frame(out)
  return(out)
}


