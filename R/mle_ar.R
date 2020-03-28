data_generator <- function(alpha, phi, sigma2, y1, n) {
  init_y_values <- function(n, acc) {
    if (n == 1) return (acc)
    else {
      vec = c(rev(acc), rep(0, length(phi) - length(acc)))
      acc = c(acc, sum(vec * phi) + alpha)
    }
    init_y_values(n-1, acc)
  }
  
  params <- c(rev(phi), alpha)
  y <- init_y_values(n=length(params), acc=c(y1))
  for (i in 1:n) {
    epsilon <- rnorm(1, 0, sqrt(sigma2))
    used_values <- c(tail(y, length(phi)), 1)
    y = c(y, sum(params * used_values) + epsilon)
  }
  return (y)
}

estimator <- function(y, modelling_params) {
  
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
  
  A = matrix(0, nrow=modelling_params+1, ncol=modelling_params+1)
  
  for (h in 0:modelling_params) {
    for (j in 0:modelling_params) {
      A[j+1, h+1] = matrix_filler(j, h)
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

  sigma2_estimator <- sum(epsilon[(modelling_params+1): length(epsilon)]^2) / (length(y) - modelling_params)

  hist(epsilon)
  
  return (list(estimated_params, mean(epsilon), sigma2_estimator))
}

#-------------------------------------------------------------


alpha <- 2
phi <-  c(0.8, .25, -0.4, 0.1) # phi_1, phi_2, phi_3. order matters!
sigma2 <- 1
modelling_params = 1

y = data_generator(alpha = alpha, phi = phi, sigma2 = sigma2, y1 = 2, n = 10000)

estimator(y=y, modelling_params=modelling_params)

parameters = list(estimator(y, 1)[1],
                  estimator(y, 2)[1],
                  estimator(y, 3)[1],
                  estimator(y, 4)[1],
                  estimator(y, 5)[1],
                  estimator(y, 6)[1],
                  estimator(y, 7)[1],
                  estimator(y, 8)[1])

#------------------------------------test
m <- 1000
result <- matrix(0, nrow=m, ncol=8)
for (i in 1:m) {
  params <- c(rev(phi), alpha)
  epsilon <- rnorm(1, 0, sqrt(sigma2))
  used_values <- c(tail(y, length(phi)), 1)
  new_value <- sum(params * used_values) + epsilon

  for (j in 1:8) {
    alpha = parameters[j][[1]][[1]][1]
    phi = parameters[j][[1]][[1]][2:(j + 1)]
    params <- c(rev(phi), alpha)
    used_values <- c(tail(y, length(phi)), 1)
    result[i, j] = abs(new_value - sum(params * used_values))
  }

}

View(result)
colMeans(result)
plot(colMeans(result), type="l")





