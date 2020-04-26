library(ggplot2)
library(MASS)
library (MCS)
library(parallel)
library(data.table)


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
resid_generator <- function(y, i_vec) {
  # y is the time series that follows an AR(p)
  # i_vec is a vector with the possible orders, e.g 1:8 if the real order is p=4
  resids = c()
  for (k in i_vec) {
    phi_k = estimator(y, k)[[1]]
    x = c(1,rev(y)[2:(k+1)])
    y_resid_k = (x %*% phi_k - rev(y)[1])^2
    resids = c(resids, y_resid_k)
  }
  
  return (resids) }

grid_search <- function(alpha, phi, sigma2, i_vec) {
  y = arima.sim(n = 2000, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
  resids = resid_generator(y, i_vec)
  phi_i = estimator(y, which.min(resids))
  prediction = c(1, rev(y)[1:(length(phi_i[[1]])-1)]) %*% phi_i[[1]]
  real = c(1, rev(y)[1:length(phi)]) %*% c(alpha, phi) + rnorm(1, 0, sigma2)
  return ((prediction - real)^2)
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
  for (k in 1:1000) 
    {
    y = arima.sim(n = 3000, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
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

