################################################################################
# Exercise 11

library(plotly)
library(pracma) # used for interp2()
source("include/MonteCarloIntegration.R")
source("include/PDEFunctions.R")

# evaluates the basis functions/eigenfunctions
v.n.m <- function(n, m, x, y, Lx=1, Ly=1){
  a <- 2/sqrt(Lx*Ly)
  b <- sin(n*pi*x/Lx)
  c <- sin(m*pi*y/Ly)
  
  return(a*b*c)
}

# evaluates the eigenvalues
beta.n.m <- function(n, m, Lx=1, Ly=1){
  return(((n/Lx)^2 + (m/Ly)^2)*pi^2)
}


N <- 8
M <- N

# x and y values on which make up the grid
x.length <- 100
y.length <- x.length
x <- seq(0, 1, length=x.length)
y <- seq(0, 1, length=y.length)

# initialise the solution variable
p.solution <- matrix(0, nrow=length(x), ncol=length(y))

# find the true solution
for(n in 1:N){
  for(m in 1:M){
    # calculate the spectral coefficients
    f.n.m <- MonteCarlo2D(user.func=function(x, y) StepFunction2D(x, y)*v.n.m(n, m, x, y))$mean
    
    # evaluate the plotting points
    for(i in 1:length(x)){
      p.solution[i, ] <- p.solution[i, ] + f.n.m*beta.n.m(n, m)^(-1)*v.n.m(n, m, x[i], y)
    }
  }
}

##### THIS IS THE BIT THAT IS PROBABLY WRONG.
# Add on noise to the true solution
p.eta <- p.solution
theta <- 1e-8*max(p.solution) # This comes from equation 48 in the notes. Should this be as in Exercise 9?

k <- 10
l <- k

for(i in 1:length(x)){
  for(j in 1:length(y)){
    p.eta[i, j] <- p.eta[i, j] + theta*v.n.m(k, l, x[i], y[j])#/2
  }
}


#################
# Solve the inverse problem
f.eta <- matrix(0, nrow=length(x), ncol=length(y))

for(n in 1:N){
  for(m in 1:M){
    # calculate the spectral coefficients
    p.n.m <- MonteCarlo2D(user.func=function(x, y){
      interp2(seq(0, 1, length=x.length), seq(0, 1, length=y.length), t(p.eta), x, y)*
      v.n.m(n, m, x, y)}, N=5000)$mean
    
    # evaluate the plotting points
    beta.n.m.evaluated <- beta.n.m(n, m, Lx=1, Ly=1)
    for(i in 1:length(x)){
      f.eta[i, ] <- f.eta[i, ] + p.n.m*beta.n.m.evaluated*v.n.m(n, m, x[i], y)
    }
  }
}


fig <- plot_ly(x=x, y=y, z=f.eta) %>% add_surface()
fig

#fields::image.plot(x, y, f.eta)





















