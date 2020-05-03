setwd("C:/Marci/CEU/ThesisDONTUSETHIS/CODE/R")
load("covidePred2.RData")
library(ggplot2)
library(plotly)
Loss
sqrt(apply(Loss, 2, sum)/nrow(Loss))
MCS # p value was 0.26

####grid search chooses m3n1 - finally some good fucking plot
params <- coeffs[[11]]
m=3
n=1
a=mn_optim(params)
plot(a$I_S + a$Q + a$R)
cases = read.csv("curve.csv")
cases$Date=ymd(cases$Date)
plot_df <- cbind(cases, "model"=a$I_S+a$Q+a$R)
p <- ggplot(data=plot_df, aes(x=Date))+
  geom_line(aes(y=model, colour="best grid search model"))+
  geom_point(aes(y=x))+
  scale_x_date(date_labels="%m-%Y")+
  geom_vline(xintercept = as.numeric(plot_df$Date[44]), 
             linetype=4, colour="black")+
  geom_hline(yintercept = 2e5, color="darkgreen")+
  geom_segment(aes(x=2020-03-01, y=2e5, xend=2020-04-04, yend=2e5), data=plot_df)
p
ggplotly(p, dynamicTicks = T)
