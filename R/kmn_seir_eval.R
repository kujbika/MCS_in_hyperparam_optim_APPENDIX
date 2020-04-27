setwd("C:/Marci/CEU/ThesisDONTUSETHIS/CODE/R")
source("kmn_seir_functions.R")
parameter_values <- c(
  k=100,
  b=0.004,
  q=1,
  N=504,
  m=2,
  sigma=0.1,
  n=3,
  gamma=1/14,
  d_I=0.8,
  tau_q = 14
)

initial_values <- c(
  S=500,
  S_Q=0,
  E1 = 1,
  E2 = 1,
  I_A1 = 1,
  I_A2 = 0,
  I_A3 = 0,
  I_S1 = 1,
  I_S2 = 0,
  I_S3 = 0,
  I_S = 1,
  I_A =1,
  I=2,
  P_I1 = 0.0,
  P_I2 = 0.0,
  P_I3 = 0.0,
  Q=0,
  R=0
)

time_values <- seq(0:100)

out = kmn_seir(parameter_values, initial_values, time_values)

write.csv(res, "C:/OwnWork/thesis/results.csv")

flu <- read.table("https://bit.ly/2vDqAYN", header = TRUE)
flu
preds <- kmn_seir(parameter_values, initial_values, flu$day)
