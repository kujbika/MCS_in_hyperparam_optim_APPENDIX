library (MCS)
data(Loss)
summary(Loss[,1:5])

rmse=function(real,...){
  forecasts=cbind(...)
  sqrt(colMeans((real-forecasts)^2))
}


library(reshape2)
library(grid)
library(gridExtra)
library(ggplot2)
library(boot)
T = 300 # = number of observations = #
n = 8 # = Number of Variables = #
rho = 0.75 # = factors AR coefficient = #
phi = 0.6 # = Variables AR coefficient = #

set.seed(1) # = seed for replication = #

# = Generate two autoregressive factors = #
f = matrix(0, T, 2)
for(i in 2:T){
  f[i, ] = rho * f[i-1, ] + rnorm(2)
}

lambda = runif(n,-1,1) # = Factor loadings = #

# = Generate variables = #
Y = matrix(0, T, n)
for(i in 2:T){
  # = Based on the first factor = #
  Y[i, 1:4] = phi * Y[i-1, 1:4] + lambda[1:4] * f[i, 1] + rnorm(4, 0, 2)
  # = Based on the second factor = #
  Y[i, 5:8] = phi * Y[i-1, 5:8] + lambda[5:8] * f[i,2] + rnorm(4, 0, 2)
}

# = Organize variable lags = #
# = y is the first (dependent) variable = #
Y = as.data.frame(embed(Y, 2)[ , c(1, 9:16)])
# = y is the first (dependent) variable = #
colnames(Y) = c("y", paste("x", 1:n, sep = ""))

# = First we will generate formulas for all the 255 equations = #
variables = colnames(Y)[-1]
combinations = lapply(1:n, function(x)combn(1:n, x))
formulas = lapply(combinations, function(y){
  apply(y, 2, function(x){
    paste("~", paste(variables[x], collapse = " + "), sep=" ")
  })
})
formulas = unlist(formulas)

# = Now we break the sample into in-sample and out-of-sample = #
in_sample = Y[1:100, ]
out_of_sample = Y[-c(1:100), ]

# = Estimate the models = #
models = lapply(formulas, function(x) lm(paste("y", x), data = in_sample))
# = Compute forecasts = #
forecasts = Reduce(cbind, lapply(models, function(x) predict(x, out_of_sample)))
colnames(forecasts) = paste("model", 1:ncol(forecasts), sep = "")

# RMSES (Mean of all models, model with all variables)
rmse(out_of_sample$y,rowMeans(forecasts),forecasts[,ncol(forecasts)])
mcs=function(Loss,R,l){
  LbarB=tsboot(Loss,colMeans,R=R,sim="fixed",l=l)$t
  Lbar=colMeans(Loss)
  zeta.Bi=t(t(LbarB)-Lbar)
  
  save.res=c()
  for(j in 1:(ncol(Loss)-1)){
    Lbardot=mean(Lbar)
    zetadot=rowMeans(zeta.Bi)
    vard=colMeans((zeta.Bi-zetadot)^2)
    
    t.stat=(Lbar-Lbardot)/sqrt(vard)
    t.max=max(t.stat)
    model.t.max=which(t.stat==t.max)
    
    t.stat.b=t(t(zeta.Bi-zetadot)/sqrt(vard))
    t.max.b=apply(t.stat.b,1,max)
    
    p=length(which(t.max<t.max.b))/R
    
    save.res=c(save.res,p)
    names(save.res)[j]=names(model.t.max)
    
    Lbar=Lbar[-model.t.max]
    zeta.Bi=zeta.Bi[,-model.t.max]
    
  }
  save.res=c(save.res,1)
  names(save.res)[j+1]=names(Lbar)
  save.p=save.res
  for(i in 2:(length(save.res)-1)){
    save.p[i]=max(save.res[i-1],save.res[i])
  }
  aux=match(colnames(Loss),names(save.p))
  save.p=save.p[aux]
  save.res=save.res[aux]
  return(list(test=save.p,individual=save.res))
}
# = generate 2 subsamples = #
out1 = out_of_sample[1:100, ]
out2 = out_of_sample[-c(1:100), ]
f1 = forecasts[1:100, ]
f2 = forecasts[-c(1:100),]

# = Compute MCS = #
set.seed(123) # = Seed for replication = #
MCS=mcs( Loss=(f1-out1$y)^2, R=1000, l=3 )$test # = R block bootstrap samples with block length = 3  = #
# = If you want to see the p-values write (print(sort(MCS)) = #
# = Selected models in the 50% MCS  = #
selected=which(MCS>0.1) # significance level

# RMSE (Mean of all models, model with all variables, MCS models)
rmse(out2$y,rowMeans(f2),f2[,ncol(f2)],rowMeans(f2[,selected]))
