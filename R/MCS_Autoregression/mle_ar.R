load("MLE_AR.RData")
library(ggplot2)
library(MASS)
library (MCS)
library(parallel)
library(data.table)
library(plotly)
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
# for example run the following
y = arima.sim(n = 100, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
estimator(y, 4)
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
# for example run the following
y = arima.sim(n = 10, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
pred.n.ahead(y, 1, estimator(y, 3)[[1]])

#-------------------------------------------------------------
resid_generator <- function(y, i_vec, n.ahead) {
  ### y is the autoregression
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
#for example run the following
y = arima.sim(n = 100, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
a = resid_generator(y, 1:8, 3)
a[[1]]
a[[2]]

#different losses and MCS
n.ahead=6
y = arima.sim(n = 1000, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
df = resid_generator(y, 1:8, n.ahead)
d=apply(df[[1]], 2, mean)
cl=makeCluster(detectCores())
MCSprocedure(df[[1]], alpha=0.1, cl=cl)
stopCluster(cl)
names(df[[1]])=paste0(1:8)
df2 <- melt(cbind(1:n.ahead, df[[1]]),  id.vars = "1:n.ahead", variable.name = 'modeled order')
p4=ggplot(df2, aes(`1:n.ahead`,value)) + geom_line(aes(colour = `modeled order`))+
  labs(x='n ahead', y="loss")
p4
ggplotly(p4, dynamicTicks=T)
#--------------------------------------
grid_search <- function(y, n.ahead, i_vec) {
  ### choses the parameter that has produced the minimal loss
  ### returns the best coefficients, and all the others as well
  df_list = resid_generator(y, i_vec, n.ahead)
  best = which.min(apply(df_list[[1]], 2, sum))
  return(list(df_list[[2]][[best]], df_list))
}
#for example run the following
y = arima.sim(n = 1000, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
a = grid_search(y, 40, 1:8)
a[[1]]
apply(a[[2]][[1]], 2, mean) #indeed, grid search chooes the minimum one


forecast_errors <- function(y, n.ahead, test.period, i_vec, real.order) {
  ### final comparison based on grid search and MCS predictions
  ### i_vec is the grid in increasing order
  u = head(y, length(y) - test.period)
  parms = grid_search(u, n.ahead, i_vec)
  forecast_grid = pred.n.ahead(u, test.period, parms[[1]][[1]])
  error_grid = (tail(forecast_grid, test.period) - tail(y, test.period))^2
  mcs_i = tail(i_vec, sum(i_vec >= real.order))
  forecast_MCS = vector(mode="list", length=length(mcs_i))
  names(forecast_MCS) = mcs_i
  for (k in mcs_i) {
    forecast_MCS[[k-real.order+1]]= pred.n.ahead(u, test.period, parms[[2]][[2]][[k-real.order+1]][[1]])
  }
  forecast_MCS = as.data.frame(forecast_MCS)
  error_MCS = tail(apply(forecast_MCS, 2, function(x) x-y), test.period)
  error_MCS = apply(error_MCS, 1, function(x) mean(x))^2
  return(list(mean(error_grid), mean(error_MCS)))
  
}
#for example run
n.ahead = 10
test.period=1
i_vec=1:8
real.order = length(phi)
y = arima.sim(n = 5000, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
forecast_errors(y, 1, 1, 1:8, real.order)


#--------------------------------
#actual research

real.order = length(phi)
results = data.frame(matrix(0, nrow=250, ncol=2))
names(results) = c("grid search loss", "MCS loss")
#mean comparison
for (l in 1:nrow(results)) {
  comparison = vector(mode="list", length=2)
  errs = as.data.frame(comparison)
  for (j in 1:1000) {
    y = arima.sim(n = 1000, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
    err = as.data.frame(forecast_errors(y, 1, 10, 1:8, real.order))
    names(err)=c("grid search loss", "MCS loss")
    errs = rbind(errs, as.data.frame(err))
  }
  results[l,]=apply(errs, 2, mean)
  if(l %% 10==0) print(l)
}
t = mean(results$MCS-results$grid)/(sqrt(var((results$MCS-results$grid))) / sqrt(250))
t
pt(t, df = 249)
t = mean(results$MCS-sigma2)/(sqrt(var((results$MCS))) / sqrt(250))
pt(t, df = 249)

p = ggplot(data = results, aes(x=`grid search loss`, y = `MCS loss`))+
  geom_density_2d()+
  geom_point()+
  geom_abline()
ggplotly(p, dynamicTicks = T)
p

#actual loss comparison
for (j in 1:10000) {
  y = arima.sim(n = 1000, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
  err = as.data.frame(forecast_errors(y, 1, 10, 1:8, real.order))
  names(err)=c("grid search loss", "MCS loss")
  errs = rbind(errs, as.data.frame(err))
}
t = mean(errs$`MCS loss`-errs$`grid search loss`)/(sqrt(var((errs$`MCS loss`-errs$`grid search loss`))) / sqrt(10000))
t
pt(t, df =9999)

errs$diff = errs$`MCS loss` < errs$`grid search loss`
pbinom(sum(errs$diff), 10000, 0.5)

df = cbind(1:10000, errs[,c(1,2)])
names(df)[1] = "Simulation number"
df <- melt(df ,  id.vars = 'Simulation number', variable.name = 'Loss')
ggplot(data= df, aes(value, after_stat(count), colour=Loss, fill=Loss)) + 
  geom_density(position='stack', alpha=0.6)+xlim(c(0,10))


ggplot(df, aes(x=`Simulation number`))+geom_density(aes(y=value, colour=Loss))
p2 = ggplot(data = errs, aes(x=`grid search loss`, y = `MCS loss`))+
  geom_density2d()+
  geom_abline()
ggplotly(p2, dynamicTicks = T)
p2
