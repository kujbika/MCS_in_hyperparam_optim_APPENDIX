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
y = arima.sim(n = 250, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
auto.arima(y, stationary=T, start.p=4, start.q=0)

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
  preds=c()
  for (k in i_vec) {
    phi_k = estimator(y, k)[[1]]
    x = c(1,rev(y)[2:(k+1)])
    y_resid_k = (x %*% phi_k - rev(y)[1])^2
    resids = c(resids, y_resid_k)
    p=c(1, rev(y)[1:(length(phi_k)-1)]) %*% phi_k
    preds=c(preds, p)
  }
  
  return (list(resids, preds)) }



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
  preds = data.frame("pred_1" = NA,
                    "pred_2" = NA,
                    "pred_3" = NA,
                    "pred_4" = NA,
                    "pred_5" = NA,
                    "pred_6" = NA,
                    "pred_7" = NA,
                    "pred_8" = NA)
  for (k in 1:1000) 
    {
    y = arima.sim(n = 1001, list(ar=c(phi)), sd=sqrt(sigma2), mean=alpha)
    
    resids = resid_generator(y[1:(length(y)-1)], 1:8)
    loss = rbind(loss, c(resids[[1]]))
    outsample = (resids[[2]] - tail(y, 1))^2
    preds = rbind(preds, c(outsample))

  }
  Loss=loss[-1,]
  preds=preds[-1,]
  grid_search_loss=c()
  for (i in 1:nrow(Loss)) {
    grid_search_loss=c(grid_search_loss, preds[i, which.min(Loss[i, ])])
  }
  
  MCS_loss=(apply(preds, 2, sum)/nrow(preds))[4:8]
  
  return (list(grid_search_loss, MCS_loss))
}

res=data.frame(matrix(0, nrow=500, ncol=2))
names(res)=c("grid search loss", "MCS loss")
for (i in 1:nrow(res)){
  a=loss_values(alpha, phi, sigma2)
  res[i, 1]=a[1]
  res[i, 2]=a[2]
}

grid_search_loss=c()
for (i in 1:nrow(Loss)) {
  grid_search_loss=c(grid_search_loss, preds[i, which.min(Loss[i, ])])
}

a = data.frame(colMeans(preds))
a$model=rownames(a)
p=ggplot(data=a, aes(y=loss, x=model, group=1))+geom_line()

mean(grid_search_loss)
a=MCSprocedure(Loss, alpha=0.15, cl=cl)
ggplot(data=preds)+

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
MCS <- MCSprocedure(Loss=a,alpha=0.1,B=5000,statistic='Tmax',cl=cl)
MCS

