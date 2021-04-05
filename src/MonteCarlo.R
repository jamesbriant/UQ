# Exercise 1
# part 1)

drawLambda <- function(n, epsilon){
  return(runif(n, 1-epsilon, 1+epsilon))
}

u <- function(t, lambda, u0=2){
  # analytical system solution
  return(u0*exp(-t*lambda))
}

m <- function(epsilon, t, u0=2){
  # analytical mean
  a <- -u0/(2*epsilon*t)
  b <- exp(-t*(1+epsilon)) - exp(-t*(1-epsilon))
  
  return(a*b)
}

t <- 5
epsilon <- 0.3
J <- 1000

lambdas <- drawLambda(J, epsilon)
solutions <- u(t, lambdas)

mean(solutions)
var(solutions)

m(epsilon, t)

# part 2)
hist(solutions)

# part 3)

t <- 5
epsilon <- 0.3
J <- 100*(1:10)
variance.estimate <- numeric(length(J))

for (j in 1:length(J)){
  lambdas <- drawLambda(J[j], epsilon)
  variance.estimate[j] <- (mean(u(t, lambdas)) - m(epsilon, t))^2
}

plot(J, log(variance.estimate))

############################################################################
# Exercise 2

library(mvtnorm)

B <- 12
t <- 5
J <- 10e5
epsilon <- 0.3

lambda.B <- seq(1-epsilon, 1+epsilon, length=B)






















