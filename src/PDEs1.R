library(plotly)

source("include/PDEFunctions.R")
source("include/MonteCarloIntegration.R")

#####################################################################
# Exercise 2

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


N <- 10
M <- N

x <- seq(0, 1, length=300)
y <- seq(0, 1, length=300)
z <- matrix(0, nrow=length(x), ncol=length(y))

for(n in 1:N){
  for(m in 1:M){
    # calculate the spectral coefficients
    f.n.m <- MonteCarlo2D(user.func=function(x, y) StepFunction2D(x, y)*v.n.m(n, m, x, y))$mean
    
    # evaluate the plotting points
    for(i in 1:length(x)){
      z[i, ] <- z[i, ] + f.n.m*beta.n.m(n, m)^(-1)*v.n.m(n, m, x[i], y)
    }
  }
}

fig <- plot_ly(x=x, y=y, z=z) %>% add_surface()
fig


##################################################
# Exercise 5

I.func <- function(x, L=1){
  if(x < L/2 && x > 0){
    return(x)
  }
  if(x < 1 && x >= L/2){
    return(L-x)
  }
  return(0)
}

t.star <- 0.03

T.x.t <- function(x, t, kappa, L=1, N=100){
  output <- numeric(length(x))
  
  for(n in 1:N){
    a <- 4*L/((n*pi)^2)
    b <- sin(n*pi/2)*sin(n*pi*x/L)
    c <- exp(-(n*pi/L)^2*kappa*t)
    output <- output + a*b*c
  }
  
  return(output)
}

kappas.N <- 100
kappas <- rlnorm(kappas.N, meanlog=log(3), 0.25)
colours.rainbow <- rainbow(kappas.N)

x <- seq(0, 1, length=200)

plot(x, T.x.t(x, t.star, kappas[1], N=10), type="l", col=colours.rainbow[i], ylim=c(0, 0.3))
for(i in 2:kappas.N){
  lines(x, T.x.t(x, t.star, kappas[i], N=10), col=colours.rainbow[i])
}

##########################################
# Exercise 7

library(fields)

GPMatern2D <- function(N, nu, tau, Lx=1, Ly=1, Nx=300, Ny=300, I=100, J=100){
  # N:    number of samples
  # nu:   matern smootheness param (0.5=exponential)
  # tau:  matern scale param
  # Nx:   discretisation points on x-axis
  # Ny:   discretisation points on y-axis
  # Lx:   x-axis domain, [0, Lx]
  # Ly:   y-axis domain, [0, Ly]
  # I:    x-axis grid, initialisation size
  # J:    y-axis grid, initialisation size
  
  # Increase I and J if circulantEmbedding() complains about "Weight function has negtive values".
  
  grid <- list(x=seq(0, I, length=100),
               y=seq(0, J, length=100))
  obj <- circulantEmbeddingSetup(grid, 
                                 cov.args=list(Covariance="Matern",
                                               smoothness=nu,
                                               theta=tau
                                 )
  )
  GRF1 <- circulantEmbedding(obj)
}

nu <- 2.5
tau <- 0.5

I <- 8
J <- I

grid <- list(x=seq(0, I, length=200),
             y=seq(0, J, length=200))
obj <- circulantEmbeddingSetup(grid, 
                               cov.args=list(Covariance="Matern",
                                             smoothness=nu,
                                             theta=tau
                                             )
                               )
GRF1 <- circulantEmbedding(obj)
image.plot(grid[[1]][1:50], grid[[2]][1:50], GRF1[1:50, 1:50])





###############################################################################

nu <- 10
tau <- 0.5

I <- 10
J <- 15

grid <- list(x=seq(0, I/J, length=200),
             y=seq(0, J/I, length=200))
obj <- circulantEmbeddingSetup(grid, 
                               cov.args=list(Covariance="Matern",
                                             smoothness=nu,
                                             range=tau/(I*J)
                                             )
                               )
GRF1 <- circulantEmbedding(obj)
image.plot(J*grid[[1]]/(I), I*grid[[2]]/(J), GRF1)






###############################################################################

set.seed(2021)
par(mfrow=c(1,2))

##########################
# 1

nu <- 0.5
tau <- 0.5

I <- 2
J <- 2

grid <- list(x=seq(0, I, length=80),
             y=seq(0, J, length=80))
obj <- circulantEmbeddingSetup(grid, 
                               cov.args=list(Covariance="Matern",
                                             smoothness=nu,
                                             range=tau
                               )
)
GRF1 <- circulantEmbedding(obj)
image.plot(grid[[1]], grid[[2]], GRF1)


##########################
# 2

set.seed(2021)
nu <- 0.5
tau <- 0.5

I <- 2
J <- 2

grid <- list(x=seq(0, I, length=100),
             y=seq(0, J, length=100))
obj <- circulantEmbeddingSetup(grid, 
                               cov.args=list(Covariance="Matern",
                                             smoothness=nu,
                                             range=tau
                               )
)
GRF1 <- circulantEmbedding(obj)
image.plot(grid[[1]], grid[[2]], GRF1)
par(mfrow=c(1,1))




#################################################################################
library(plgp)





# kernel function
rbf_D <- function(X,l=1, eps = sqrt(.Machine$double.eps) ){
  D <- plgp::distance(X)
  Sigma <- exp(-D/l)^2# + diag(eps, nrow(X))
}
# number of samples
nx <- 30
x <- seq(0,1,length=nx)
#grid of pairwise values
X <- expand.grid(x, x)
# compute squared exponential kernel on pairwise values
Sigma <- rbf_D(X,l=2)

# sample from multivariate normal with mean zero, sigma = sigma
Y <- MASS::mvrnorm(1,rep(0,dim(Sigma)[1]), Sigma)

# plot results
pp <- data.frame(y=Y,x1=X[,1],x2=X[,2])
z <- matrix(pp$y, ncol=30, nrow=30, byrow=FALSE)

fig <- plot_ly(x=x, y=x, z=z) %>% add_surface()
fig




















############################################################################
# Karhunen-Loeve attempt

source("include/GPFunctions2D.R")
library(plotly)

I <- 50
J <- 50
M <- I*J

x <- seq(0, 1, length=I)
y <- seq(0, 1, length=J)

tic2 <- Sys.time()
z <- GenerateGP2D(1, M=M, I=I)
Sys.time() - tic2

tic2 <- Sys.time()
z <- GenerateGP2D(4, M=M, I=I)
Sys.time() - tic2


par(mfrow=c(2,2))
persp(x, y, z[[1]])
persp(x, y, z[[2]])
persp(x, y, z[[3]])
persp(x, y, z[[4]])
par(mfrow=c(1,1))


fig <- plot_ly(x=x, y=y, z=z[[1]]) %>% add_surface()
fig








