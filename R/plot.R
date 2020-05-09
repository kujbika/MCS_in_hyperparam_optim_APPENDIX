library(lubridate)
library(data.table)
library(ggplot2)
cases = read.csv("curve.csv")
cases$Date=ymd(cases$Date)
names(cases)=c("Date", "Total confirmed")
df <- melt(cases ,  id.vars = "Date", variable.name = 'series')

# plot on same grid, each series colored differently -- 
# good if the series have same scale
p2=ggplot(df, aes(Date,value)) + geom_point(aes(colour = series))+ 
  scale_x_date(date_labels="%Y-%m-%d")+
  scale_y_continuous(name="Covid-19 cases, New York state, U.S.A")
p2
ggplotly(p2, dynamicTicks = T)
