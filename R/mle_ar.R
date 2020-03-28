alpha <- 3
phi <-  c(0.8, .25, -0.4) # phi_1, phi_2, phi_3. order matters!
modelling_params = 3
params <- c(rev(phi), alpha)
sigma2 <- 1
y1 = 2
init_y_values <- function(n=length(params), acc=c(y1)) {
  if (n == 1) return (acc)
  else {
    vec = c(rev(acc), rep(0, length(phi) - length(acc)))
    acc = c(acc, sum(vec * phi) + alpha)
  }
  init_y_values(n-1, acc)
}


y <- init_y_values()
  
n = 20000+length(params)
 
for (i in 1:n) {
  epsilon <- rnorm(1, 0, sqrt(sigma2))
  used_values <- c(tail(y, length(phi)), 1)
  y = c(y, sum(params * used_values) + epsilon)
}

y = y[floor(n/2): n]
plot(y, type='l')



matrix_filler <- function(j, h) {
  # h and j goes from zero to i
  # j stands for rows, h stands for columns
  sum1 <- 0
  sum2 <- 0
  if (h == 0){
    if (j == 0){
      return (1)
    } else {
      for (k in (modelling_params+1):length(y)){
        sum1 = sum1 + y[k-j]
        sum2 = sum2 + y[k-j]*y[k-j]
      }
      return (sum1/sum2)
    }
  } else if (j == 0){
    for (k in (modelling_params+1):length(y)) {
      sum1 = sum1 + y[k-h]
    }
    return (sum1/(length(y) - modelling_params-1)) # LEHET HOGY I+1ET KELL KIVONNI
  } else {
    for (k in (modelling_params+1):length(y)) {
      sum1 = sum1 + y[k-h] * y[k-j]
      sum2 = sum2 + y[k-j] * y[k-j]
    }
    return (sum1/sum2)
  }
}

A = matrix(0, nrow=modelling_params+1, ncol=modelling_params+1)

for (h in 0:modelling_params) {
  for (j in 0:modelling_params) {
    A[j+1, h+1] = matrix_filler(j, h)
  }
}

vector_filler <- function(j) {
  s <- 0
  squared_sum <- 0
  if (j == 0){
    return (sum(y[(modelling_params+1): length(y)]) / (length(y) - modelling_params-1)) # LEHET HOGY I+1ET KELL KIVONNI
  } else {
    for (k in (modelling_params+1):length(y)) {
      s = s + y[k] * y[k-j]
      squared_sum = squared_sum + y[k-j] * y[k-j]
    }
    return (s/squared_sum)
  }

}
b = rep(0, modelling_params+1)
for (j in 0:modelling_params) {
  b[j+1] = vector_filler(j)
}

estimated_params <- solve(A) %*% b


# now, lets create the error terms
epsilon <- rep(0, length(y))
for (d in length(y): (modelling_params+1)) {
  used_values = c(1, rev(y[(d-modelling_params):(d-1)]), y[d])
  epsilon[d] = used_values[modelling_params+2] - used_values[1:(modelling_params+1)] %*% estimated_params
}
mean(epsilon)
var(epsilon)

sigma2_estimator <- sum(epsilon[(modelling_params+1): length(epsilon)]^2) / (length(y) - modelling_params)
estimated_params
sigma2_estimator

hist(epsilon)
