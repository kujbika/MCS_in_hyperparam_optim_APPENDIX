load("covidPred_secondoptim.RData")
library(ggplot2)
library(plotly)
library(lubridate)
#the Loss and MCS objects are in the RData file
Loss
min(sqrt(apply(Loss, 2, sum)/nrow(Loss)))
mean((apply(Loss, 2, sum)/nrow(Loss))[c("m2n1", "m3n2", "m3n3", "m3n4", "m3n5",
                                    "m5n2", "m5n3", "m5n5")])
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
  scale_y_continuous(name="Total confirmed cases, N.Y state, U.S.A")+
  geom_vline(xintercept = as.numeric(cases$Date[44]), 
             linetype=4, colour="black")+
  geom_vline(xintercept = as.numeric(cases$Date[59]), 
             linetype=4, colour="black")
 
p

ggplotly(p, dynamicTicks = T)


cases2= as.data.frame(lapply(cases[,-1], diff, lag=1))
cases2 = cbind(cases[-1,]$Date, cases2)
names(cases2)[1]="Date"
p2 <- ggplot(data=cases2, aes(x=Date))+
  geom_line(aes(y=m3n4, colour="Grid search"))+
  geom_line(aes(y=MCS.prediction, colour="MCS0.90"))+
  geom_line(aes(y=Actual, colour="Actual"))+
  scale_colour_manual("",
                      breaks=c("Grid search", "MCS0.90","Actual"),
                      values=c("red", "blue", "black"))+
  scale_x_date(date_labels="%Y-%m-%d")+
  scale_y_continuous(name="New confirmed cases, N.Y state, U.S.A")+
  geom_vline(xintercept = as.numeric(cases$Date[44]), 
             linetype=4, colour="black")+
  geom_vline(xintercept = as.numeric(cases$Date[59]), 
             linetype=4, colour="black")

p2

ggplotly(p2, dynamicTicks = T)

#Some explanatory figures about the best fitting model
best_model = c(3,4)
source(path(best_model[1],best_model[2]))
parms <- coeffs[[ 5*(best_model[1]-1) + best_model[2] ]]
a=mn_dataframe(parms)
cases = read.csv("curve.csv")
cases$Date=ymd(cases$Date)
names(cases) = c("Date", "Actual")
cases = cbind(cases, a)
p_S <- ggplot(data=cases, aes(x=Date))+guides(color = FALSE, size = FALSE)+
  geom_line(aes(y=S, colour="Susceptible"))+
  scale_x_date(date_labels="%Y-%m-%d")+
  scale_y_continuous(name="Susceptible class, N.Y state, U.S.A")
p_S
ggplotly(p_S, dynamicTicks = T)
p_S_Q <- ggplot(data=cases, aes(x=Date))+guides(color = FALSE, size = FALSE)+
  geom_line(aes(y=S_Q, colour="Susceptible quarantine"))+
  scale_x_date(date_labels="%Y-%m-%d")+
  scale_y_continuous(name="Susceptible quarantine class, N.Y state, U.S.A")
p_S_Q
ggplotly(p_S_Q, dynamicTicks = T)

p_I <- ggplot(data=cases, aes(x=Date))+
  geom_line(aes(y=I_A, colour="Asymptotic infectious"))+
  geom_line(aes(y=I_S, colour="Symptotic infectious"))+
  scale_x_date(date_labels="%Y-%m-%d")+
  scale_y_continuous(name="Infectious class, N.Y state, U.S.A")
p_I
ggplotly(p_I, dynamicTicks = T)

p_Q <- ggplot(data=cases, aes(x=Date))+guides(color = FALSE, size = FALSE)+
  geom_line(aes(y=Q, colour="Asymptotic infectious"))+
  scale_x_date(date_labels="%Y-%m-%d")+
  scale_y_continuous(name="Infectious quarantine class, N.Y state, U.S.A")
p_Q
ggplotly(p_Q, dynamicTicks = T)

p_E <- ggplot(data=cases, aes(x=Date))+guides(color = FALSE, size = FALSE)+
  geom_line(aes(y=E1+E2+E3, colour="Asymptotic infectious"))+
  scale_x_date(date_labels="%Y-%m-%d")+
  scale_y_continuous(name="Exposed class, N.Y state, U.S.A")
p_E
ggplotly(p_E, dynamicTicks = T)


