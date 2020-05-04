setwd("C:/Marci/CEU/ThesisDONTUSETHIS/CODE/R")
load("covidPred_secondoptim.RData")
library(ggplot2)
library(plotly)
library(lubridate)
Loss
min(sqrt(apply(Loss, 2, sum)/nrow(Loss)))
MCS # p value was 0.15

####grid search chooses and mcs preds
hyp_parms <- list("m2n1"=c(2,1), "m3n2"=c(3,2), "m3n3"=c(3,3),
                  "m3n4"=c(3,4), "m3n5"=c(3,5), "m5n2"=c(5,2),
                  "m5n3"=c(5,3), "m5n5"=c(5,5))

cases = read.csv("curve.csv")
cases$Date=ymd(cases$Date)
for (i in hyp_parms) {
  source(path(i[1],i[2]))
  parms <- coeffs[[ 5*(i[1]-1) + i[2] ]]
  a=mn_dataframe(parms)
  cases=cbind(cases, a$I_S+a$Q+a$R)
}
names(cases)=c("Date", "Actual", names(hyp_parms))
cases$`MCS prediction`=apply(cases[,c(3:ncol(cases))], 1, sum)/(ncol(cases)-2)
summary(cases)



p <- ggplot(data=cases, aes(x=Date))+
  geom_line(aes(y=m3n4, colour="Grid search"))+
  geom_line(aes(y=`MCS prediction`, colour="MCS0.90"))+
  geom_point(aes(y=Actual, colour="Actual"))+
  scale_colour_manual("",
                      breaks=c("Grid search", "MCS0.90","Actual"),
                      values=c("red", "blue", "black"))+
  scale_x_date(date_labels="%Y-%m-%d")+
  scale_y_continuous(name="Epidemic curve, N.Y state, U.S.A")+
  geom_vline(xintercept = as.numeric(cases$Date[44]), 
             linetype=4, colour="black")+
  geom_vline(xintercept = as.numeric(cases$Date[59]), 
             linetype=4, colour="black")
 
p

ggplotly(p, dynamicTicks = T)

matplot(t(cases), type = 'l', col = seq_len(nrow(df1)), xaxt = "n")
legend("top", row.names(df1), col = seq_len(nrow(df1)), fill = seq_len(nrow(df1)))
axis(1, at = seq_along(df1), labels = colnames(df1))
