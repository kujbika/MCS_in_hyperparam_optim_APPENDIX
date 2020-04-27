require("deSolve")
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


prob_gen <- function(m, sigma, E_m_lagged, tau_d, n, gamma, i) {
  return(m * sigma * E_m_lagged * exp(-n*gamma*tau_d)*((n*gamma*tau_d)^(i-1))/(i-1))
}

kmn_seir <- function(parameter_values, initial_values, time_values) {
  out = dede(
    y=initial_values,
    times=time_values,
    func=kmn_seir_equations,
    parms = parameter_values
  )
  return(as.data.frame(out))
}


out = kmn_seir()





