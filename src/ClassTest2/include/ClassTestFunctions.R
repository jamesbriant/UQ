u <- function(lambda){
  beta.m <- c(20, 20, 5, 5)
  lambda.m <- c(-0.7, -0.3, 0.3, 0.7)
  
  sum.m <- 0
  
  for(m in 1:4){
    numerator <- (-1)^(m+1)
    denominator <- 1 + exp(-2*beta.m[m]*(lambda - lambda.m[m]))
    
    sum.m <- sum.m + numerator/denominator
  }
  
  return(sum.m)
}

dwigner <- function(x, R=1){
  a <- 2/(pi*R^2)
  b <- sqrt(R^2 - x^2)
  
  return(a*b)
}

rwigner <- function(n, R=1){
  Y <- rbeta(n, 1.5, 1.5)
  return(R*(2*Y - 1))
}