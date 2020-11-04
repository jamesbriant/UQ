library(deSolve)
library(mvtnorm)

n <- 10
alphas <- seq(5, 10, length=n) # range of alpha values


# initial conditions
x0 <- 5
y0 <- 8
intial.conditions <- c(x0, y0)

# set up time interval with many points for more accurate results
tt <- (1:2000)/100

# build the LV equations
odes <- function(t, data, params){
  x <- data[1]
  y <- data[2]
  
  alpha <- params[1]
  
  xprime <- alpha*x - 0.5*x*y
  yprime <- 0.6*x*y - y
  
  return(list(c(xprime, yprime)))
}

# prep vector for solutions at x(20)
x20.estimates <- numeric(n)

# solve equations for each alpha, saving x(20) each time
for(i in 1:length(alphas)){
  solution <- as.data.frame(ode(intial.conditions, tt, odes, c(alphas[i])))
  names(solution) <- c("time", "x", "y")
  x20.estimates[i] <- solution$x[length(tt)]
}

# create the data frame using the solution results
D <- data.frame(alphas, x20.estimates)
colnames(D) <- c("alpha", "zeta")

# function for calculating the variance
calcSigma <- function(X1, X2){
  # Function takes in a discretised grid
  #input X1 x X2 and computes covariance matrix using covariance function k
  
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)){
    for (j in 1:ncol(Sigma)){
      Sigma[i,j] <- exp(-0.5*(abs(X1[i] - X2[j]))^2) + 0.15*(i==j)
    }
  }
  
  return(Sigma)
}

########## RUN THIS BLOCK BELOW FOR THE FIRST PLOT ####################

# calculate the variance matrices required for plotting
alpha.plot <- seq(5, 10, length=50)
k.aa <- calcSigma(D$alpha, D$alpha)
k.aap <- calcSigma(D$alpha, alpha.plot)
k.apa <- calcSigma(alpha.plot, D$alpha)
k.apap <- calcSigma(alpha.plot, alpha.plot)

# calculate the best solution
m.star <- k.apa %*% solve(k.aa + 0.15*diag(1, ncol(k.aa))) %*% D$zeta
k.star <- k.apap - k.apa %*% solve(k.aa + 0.15*diag(1, ncol(k.aa))) %*% k.aap

# plot the solution
plot(alpha.plot, m.star, xlab="alpha", ylab="Prey Population (zeta)", main="Prey Population Estimates at t=20 for different alpha", ylim=c(-1.5, 7))
lines(alpha.plot, m.star)
points(D$alpha, D$zeta, pch=16, col="red")

# calculate and plot the confidence bands
y.up <- m.star + 2.5758*diag(k.star)^0.5
y.low <- m.star - 2.5758*diag(k.star)^0.5
lines(alpha.plot, y.low, col="dark grey", lty=2)
lines(alpha.plot, y.up, col="dark grey", lty=2)
# ignore the following warnings
legend(8, 6, legend=c("best estimate", "99% confidence band", "data points"), col=c("black", "grey", "red"), lty=c(1, 2, 0), pch=c(26, 26, 16))

####### END OF FIRST BLOCK ###########



########## RUN THIS BLOCK BELOW FOR THE SECOND PLOT ####################

# Now calculate the distribution of g(alpha=7.5)

# calculate the variance matrices required
alpha.plot <- 7.5
k.aa <- calcSigma(D$alpha, D$alpha)
k.aap <- calcSigma(D$alpha, alpha.plot)
k.apa <- calcSigma(alpha.plot, D$alpha)
k.apap <- calcSigma(alpha.plot, alpha.plot)

# calculate the best solution
m.star <- k.apa %*% solve(k.aa + 0.15*diag(1, ncol(k.aa))) %*% D$zeta
k.star <- k.apap - k.apa %*% solve(k.aa + 0.15*diag(1, ncol(k.aa))) %*% k.aap

hist(rnorm(1000, m.star, k.star), breaks=seq(1, 3.1, 0.1), xlab="Prey Population - g(alpha=7.5)", main="Distribution of Prey Population drawn from 1000 samples of g(7.5)=N(2.10, 0.234)")

####### END OF SECOND BLOCK ###########


# The first plot shows the best estimate of the prey population at time 20 for different values of alpha on [5,10].
# The population is estimated to be largest at t=20 when alpha is approximately 6.8, however the error here is large.

# The histogram shows the distribution of the predicted population when alpha=7.5 from 1000 samples. 
# The population mean is approximate 2.1 and the variance is approximately 0.234.
# This is repre









