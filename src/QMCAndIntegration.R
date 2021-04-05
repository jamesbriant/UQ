source("include/QMC.R")
source("include/ExponentialExample.R")
source("include/Quadrature.R")

############################################################################

J <- 1000
VDC.sequence <- VDCSequence(J, base=2)

lambda.dagger <- 2
epsilon       <- 0.5
lambda.sample <- lambda.dagger - epsilon + 2*epsilon*VDC.sequence

t <- 5
QMC.sample <- u(t, lambda.sample)
mean(QMC.sample)

#################
# Find the order of convergence

## Use of MSE makes no sense as the values are not random, they're deterministic

Js <- seq(20, 50000, by=20)
I <- length(Js)
square.error <- numeric(I)
absolute.error <- numeric(I)
exact.sol <- m(t, lambda.dagger, epsilon)

pb <- txtProgressBar(min = 0, max = I, style = 3)
for(i in 1:I){
  lambda.sample <- lambda.dagger - epsilon + 2*epsilon*VDCSequence(Js[i], base=2)
  square.error[i] <- (exact.sol - mean(u(t, lambda.sample)))^2
  absolute.error[i] <- abs(exact.sol - mean(u(t, lambda.sample)))
  setTxtProgressBar(pb, i)
}

plot(Js, absolute.error, type="l", log="xy")
lines(Js, log(Js)/(Js*10000), col="blue")






#############################################################################
# Exercise 4

f <- function(x) u(t, x)

SR <- SimpsonRule(f, lambda.dagger - epsilon, lambda.dagger + epsilon)
SR
CSR <- CompositeSimpsonRule(f, 5, lambda.dagger - epsilon, lambda.dagger + epsilon)
CSR

############################################################################
# Exercise 5

x <- 2:21
estimates <- sapply(x, GaussQuadrature, user.func=f, lambda.dagger=lambda.dagger, epsilon=epsilon)
estimates
plot(x, estimates, type="l")
abline(h=SR, col="red")
abline(h=CSR, col="blue")


###############################################################################
# Exercise 6

g <- function(x, beta=1){
  
}