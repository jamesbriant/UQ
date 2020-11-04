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

# print out the estimates at x(20)
x20.estimates

# create the data frame using the solution results
D <- data.frame(alphas, x20.estimates)
colnames(D) <- c("alpha", "zeta")






























