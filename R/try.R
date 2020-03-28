library(MCS)
library(rugarch)
data("STOXXIndexesRet")
ret <- STOXXIndexesRet[, "SXA1E"]

models <- c('sGARCH', 'eGARCH', 'gjrGARCH', 'apGARCH', 'csGARCH')
distributions <- c("norm", "std", "ged", "snorm", "sstd", "sged", "jsu", "ghyp")
spec.comp <- list()
for(m in models){
  for(d in distributions){
    spec.comp[[paste(m, d, sep = "-")]] <- ugarchspec(
      mean.model = list(armaOrder = c(0,0)), 
      variance.model = list(model = m, garchOrder =c(1,1)),
      distribution.model = d)
  }
}
specifications <- names(spec.comp)

roll.comp <- list()
for(s in specifications){
  roll.comp[[s]] <- ugarchroll(spec = spec.comp[[s]], data = ret, 
                               forecast.length = 2000, refit.every = 200)
}

VaR.comp = list()
for(s in specifications){
  VaR.comp[[s]] = as.data.frame(roll.comp[[s]], which = "VaR")[, 1]
}

Loss <- do.call(cbind, lapply(specifications, 
                              function(s) LossVaR(tau = 0.01, realized = tail(ret, 2000)/100,
                                                  evaluated = VaR.comp[[s]]/100)))
colnames(Loss) <- specifications
# alpha is the significance level (1 - alpha is the confidence level)
SSM <- MCSprocedure(Loss = Loss, alpha = 0.05, B = 5000, cl = NULL, statistic = "Tmax")
