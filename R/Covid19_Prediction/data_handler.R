data <- read.csv("time_series_covid19_confirmed_US.csv")
curve <- apply(data, 2, sum)
write.csv(curve, 'curve.csv')
