#########################################################################

library(mvtnorm)

# Question 1

radialBasis <- function(r, s, t){
  a <- (s-t)^2
  b <- 2*r^2
  
  return(exp(-a/b))
}

brownianMotion <- function(s, t){
  return(min(s,t))
}

#whiteNoise <- function(){}

matern <- function(s, t){
  a <- s-t
  
  return((1+abs(a))*exp(t-s))
}

n <- 300
t <- seq(from=0, to=5, length.out = n)
m <- rep(0, n)

K <- matrix(0, n, n)
for(i in 1:n){
  K[i, ] <- sapply(0.25, radialBasis, t[i], t)
}

plot(t, rmvnorm(1, m, K), type="l")
for(i in 1:14){
  lines(t, rmvnorm(1, m, K), type="l")  
}



############-----------------------------------------------------------
# Exercise 2

###### Part 1

phi <- function(j, x){
  return(sqrt(2)*sin(j - pi/(2*x)))  
}


M <- 200
c <- rnorm(M, 0, 1)
x <- seq(from=0.001, to=1, length.out = M)

fx <- c*mapply(phi, 1:M, x)

plot(x, fx, type="l")


######## Part 2

phi <- function(j, x){
  return(sin(j*pi*x))  
}

M <- 200
c <- rnorm(M, 0, 1)
x <- seq(from=0.001, to=1, length.out = M)

fx <- c*mapply(phi, 1:M, x)

plot(x, fx, type="l")











