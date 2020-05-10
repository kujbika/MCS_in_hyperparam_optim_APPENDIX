library(ggplot2)
library(MASS)
library (MCS)
library(parallel)
library(data.table)

#parameter declarations for the AR process
#alpha is the intercept, phi is the coeff vector. Phi can be any long
#sigma2 is the variance
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

estimator <- function(y, modelling_params) {
  ###parameter estimation function based on Conditional maximum likelihoog maximization
  ###for details see https://github.com/kujbika/MCS_in_hyperparam_optim_APPENDIX/tree/master/R/MCS_Autoregression
  row_number = length(y) - modelling_params + 1
  vec = matrix(rev(y)[1 : modelling_params + 1], ncol = modelling_params)
  for (i in 3:row_number) {
    vec = rbind(vec, matrix(rev(y)[i : (i + modelling_params - 1)], ncol=modelling_params))
  }

  X = cbind(1, vec)
  X_Mp_inverse <- solve(t(X) %*% X) %*% t(X)
  phi = X_Mp_inverse %*% head(rev(y), row_number-1)
  sigma2_est =  (1/row_number) * t(head(rev(y), row_number-1) - X %*% phi) %*% (head(rev(y), row_number-1) - X %*% phi)
  return (list(phi, sigma2_est))
}
#-------------------------------------------------------------
pred.n.ahead <- function(y, n.ahead, parms) {
  ### n ahead prediction based on given parameters
  if (n.ahead == 0) {
    return(y)
  }else {
    prediction = c(1, rev(y)[1:(length(parms)-1)]) %*% parms
    return(pred.n.ahead(c(y, prediction), n.ahead - 1, parms))
  }
  
}

resid_generator <- function(y, i_vec, n.ahead) {
  ### y is the time series that follows an AR(p)
  ### i_vec is a vector with the possible orders, e.g 1:8 if the real order is p=4
  ### returns n.ahead prediction residual terms
  pars <- vector(mode="list", length=length(i_vec))
  names(pars) = i_vec
  coeffs <- vector(mode="list", length=length(i_vec))
  names(coeffs) = i_vec
  u = head(y, length(y) - n.ahead)
  for (k in i_vec) {
    parms <- estimator(u, k)
    forecast = pred.n.ahead(u, n.ahead, parms[[1]])
    coeffs[[k]] = parms
    pars[[k]] = (tail(forecast, n.ahead) - tail(y, n.ahead))^2
    }
  return (list(as.data.frame(pars), coeffs)) }

#research
n.ahead=20
y = arima.sim(n = 1000, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
df = resid_generator(y, 1:8, n.ahead)
d=apply(df[[1]], 2, mean)
cl=makeCluster(detectCores())
MCSprocedure(df[[1]], alpha=0.1, cl=cl)
stopCluster(cl)
df2 <- melt(cbind(1:n.ahead, df[[1]]),  id.vars = "1:n.ahead", variable.name = 'series')
ggplot(df2, aes(`1:n.ahead`,value)) + geom_line(aes(colour = series))

#--------------------------------------
grid_search <- function(y, n.ahead, i_vec) {
  ### choses the parameter that has produced the minimal loss
  ### returns the best coefficients, and all the others as well
  df_list = resid_generator(y, i_vec, n.ahead)
  best = which.min(apply(df_list[[1]], 2, sum))
  return(list(df_list[[2]][[best]], df_list[[2]]))
}

forecast_errors <- function(y, n.ahead, test.period, i_vec) {
  ### final comparison based on grid search and MCS predictions
  u = head(y, length(y) - test.period)
  parms = grid_search(u, n.ahead, i_vec)
  forecast_grid = pred.n.ahead(u, test.period, parms[[1]][[1]])
  forecast_MCS = vector(mode="list", length=length(i_vec))
  names(forecast_MCS) = i_vec
  for (k in i_vec) {
    forecast_MCS[[k]] =pred.n.ahead(u, test.period, parms[[2]][[k]][[1]])
  }
  forecast_MCS = as.data.frame(forecast_MCS)
  error_grid = (tail(forecast_grid, test.period) - tail(y, test.period))^2
  error_MCS = tail(apply(forecast_MCS, 2, function(x) x - y), test.period)
  error_MCS = apply(error_MCS, 1, function(x) mean(x))^2
  return(list(mean(error_grid), mean(error_MCS)))
  
}
comparison = vector(mode="list", length=2)
errs = as.data.frame(comparison)

  i=1.2
  for (j in 1:1000) {
    y = arima.sim(n = 100, list(ar=c(phi)), sd=sqrt(i), mean=alpha)
    err = as.data.frame(forecast_errors(y, 20, 16, 1:8))
    names(err)=c("grid", "MCS")
    errs = rbind(errs, as.data.frame(err))
    if(j %% 10==0) print(j)
  }
t = mean(errs$MCS-errs$grid)/sqrt(var((errs$MCS-errs$grid)))
errs$diff = errs$MCS < errs$grid
sum(errs$diff)

mins = apply(a, 1, min)
mins = data.table("minimal observed losses"= mins)
mins$`minimal observed losses`=as.numeric(mins$`minimal observed losses`)
ggplot(mins, aes(x=`minimal observed losses`)) + 
  geom_histogram(aes(y=..density..), color="darkblue", fill="lightblue")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(aes(xintercept=mean(`minimal observed losses`)),
             color="blue", linetype="dashed", size=1)
fit1 <- fitdistr(mins$`minimal observed losses`, "exponential")
ks.test(mins$`minimal observed losses`, "pexp", fit1$estimate)
fit1


