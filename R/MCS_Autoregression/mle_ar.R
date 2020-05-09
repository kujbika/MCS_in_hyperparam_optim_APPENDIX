library(ggplot2)
library(MASS)
library (MCS)
library(parallel)
library(data.table)
load("MLE_AR.RData")


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
y = arima.sim(n = 1000, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)


estimator <- function(y, modelling_params) {
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
estimator(y, 4) 
#-------------------------------------------------------------
pred.n.ahead <- function(y, n.ahead, parms) {
  
  if (n.ahead == 0) {
    return(y)
  }else {
    prediction = c(1, rev(y)[1:(length(parms)-1)]) %*% parms
    return(pred.n.ahead(c(y, prediction), n.ahead - 1, parms))
  }
  
}

resid_generator <- function(y, i_vec, n.ahead) {
  # y is the time series that follows an AR(p)
  # i_vec is a vector with the possible orders, e.g 1:8 if the real order is p=4
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

y = arima.sim(n = 1000, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
df = resid_generator(y, 1:8, 16)
cl=makeCluster(detectCores())
MCSprocedure(df[[1]], alpha=0.05, cl=cl)
stopCluster(cl)
df2 <- melt(cbind(1:16, df[[1]]),  id.vars = '1:16', variable.name = 'series')
ggplot(df2, aes(`1:16`,value)) + geom_line(aes(colour = series))


grid_search <- function(y, n.ahead, i_vec) {
  df_list = resid_generator(y, i_vec, n.ahead)
  best = which.min(apply(df_list[[1]], 2, sum))
  return(df_list[[2]][[best]])
}

grid_search_forecast_error <- function(y, n.ahead, test.period, i_vec) {
  u = head(y, length(y) - test.period)
  parms = grid_search(u, n.ahead, i_vec)
  forecast = pred.n.ahead(head(u, length(u) - n.ahead), n.ahead+test.period, parms[[1]])
  return((tail(forecast, test.period) - tail(y, test.period))^2)
  
}


loss_values <- function(alpha, phi, sigma2) {
  # alpha is the intercept of the AR(p)
  # phi is the real coeff. vector of the same AR(p)
  # this is only run for experimantla purposes, so i_vec here is alwayzs 1:8
  
  loss = data.frame("loss_1" = NA,
                    "loss_2" = NA,
                    "loss_3" = NA,
                    "loss_4" = NA,
                    "loss_5" = NA,
                    "loss_6" = NA,
                    "loss_7" = NA,
                    "loss_8" = NA)
  for (k in 1:500) 
    {
    y = arima.sim(n = 500, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
    resids = resid_generator(y, 1:8)
    loss = rbind(loss, c(resids))
  }
  return (loss[-1,])
}

a = loss_values(alpha, phi, sigma2)
cor(a)
summary(a)

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


#######mcs
number_cores <- detectCores()
cl <- makeCluster(number_cores)
MCS <- MCSprocedure(Loss=a,alpha=0.2,B=5000,statistic='Tmax',cl=cl)
MCS

