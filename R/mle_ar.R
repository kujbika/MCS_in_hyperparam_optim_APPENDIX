alpha <- 5
phi <-  c(0.8, .25, -0.4, 0.1) # phi_1, phi_2, phi_3. order matters!
isStationer <- function(phi) {
  if (min(abs(polyroot(c(1, -phi)))) <=1) {
    print('Phi does not describe a stationary process')
    phi = as.double(strsplit(readline("Define a new phi vector, with the following format: phi1,phi2,..."), ",")[[1]])
    return (isStationer(phi))
  } else return(phi)
}
phi = isStationer(phi)
sigma2 <- 1.2
y = arima.sim(n = 10000, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)

modelling_params = 4 #the order of the model, like it was unknown before

estimator <- function(y, modelling_params) {
  row_number = length(y) - modelling_params + 1
  vec = matrix(rev(y)[1 : modelling_params+1], ncol = modelling_params)
  for (i in 3:row_number) {
    vec = rbind(vec, matrix(rev(y)[i : (i + modelling_params - 1)], ncol=modelling_params))
  }

  X = cbind(1, vec)
  X_Mp_inverse <- solve(t(X) %*% X) %*% t(X)
  phi = X_Mp_inverse %*% head(rev(y), row_number-1)
  sigma2_est =  (1/row_number) * t(head(rev(y), row_number-1) - X %*% phi) %*% (head(rev(y), row_number-1) - X %*% phi)
  return (list(phi, sigma2_est))
}
estimator(y, 4) 
#-------------------------------------------------------------





