library(plotly)

source("include/PDEFunctions.R")
source("include/MonteCarloIntegration.R")

###############################################################################
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


###############################################################################
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

###############################################################################
# Exercise 7
# Karhunen-Loeve implementation

source("include/GPFunctions2D.R")
#library(plotly)

I <- 35
J <- I
N <- 4 # number of samples

x <- seq(0, 1, length=I)
y <- seq(0, 1, length=J)
z <- GenerateGP2D(N, I=I, J=J)

par(mfrow=c(2,2))
fields::image.plot(x, y, z[[1]])
fields::image.plot(x, y, z[[2]])
fields::image.plot(x, y, z[[3]])
fields::image.plot(x, y, z[[4]])
par(mfrow=c(1,1))


#fig <- plot_ly(x=x, y=y, z=z[[1]]) %>% add_surface()
#fig

########################



source("include/GPFunctions2D.R")
#library(plotly)

I <- 35
J <- I
N <- 4 # number of samples

x <- seq(0, 1, length=I)
y <- seq(0, 1, length=J)

KLDecomp1 <- GetKLDecomposition(I=I, J=J, nu=0.1, sigma2=1, tau=1) # THIS IS VERY SLOW
KLDecomp2 <- GetKLDecomposition(I=I, J=J, nu=1, sigma2=1, tau=1) # THIS IS VERY SLOW
KLDecomp3 <- GetKLDecomposition(I=I, J=J, nu=1.5, sigma2=1, tau=1) # THIS IS VERY SLOW
KLDecomp4 <- GetKLDecomposition(I=I, J=J, nu=3, sigma2=1, tau=1) # THIS IS VERY SLOW
KLDecomp5 <- GetKLDecomposition(I=I, J=J, nu=0.5, sigma2=10, tau=1) # THIS IS VERY SLOW

#samples <- GenerateSamples(N, KLDecomp)

par(mfrow=c(2,2))
fields::image.plot(x, y, GenerateSamples(1, KLDecomp1, I=I, J=J)[[1]])
fields::image.plot(x, y, GenerateSamples(1, KLDecomp2, I=I, J=J)[[1]])
fields::image.plot(x, y, GenerateSamples(1, KLDecomp3, I=I, J=J)[[1]])
fields::image.plot(x, y, GenerateSamples(1, KLDecomp4, I=I, J=J)[[1]])
par(mfrow=c(1,1))

par(mfrow=c(1,2))
contour(x, y, GenerateSamples(1, KLDecomp1, I=I, J=J)[[1]])
contour(x, y, GenerateSamples(1, KLDecomp4, I=I, J=J)[[1]])
par(mfrow=c(1,1))

#fig <- plot_ly(x=x, y=y, z=GenerateSamples(1, KLDecomp1, I=I, J=J)[[1]]) %>% add_surface()
#fig

###############

I <- 35
J <- I
N <- 1 # number of samples

x <- seq(0, 1, length=I)
y <- seq(0, 1, length=J)

nu <- c(0.5, 1.5, 2.5, 3.5)

samples <- list()
for(i in 1:4){
  KLDecomp <- GetKLDecomposition(I=I, J=J, nu=nu[i]) # THIS IS VERY SLOW

  samples[[i]] <- GenerateSamples(N, KLDecomp, I=I, J=J)[[1]]
}

par(mfrow=c(2,2))
fields::image.plot(x, y, samples[[1]])
fields::image.plot(x, y, samples[[2]])
fields::image.plot(x, y, samples[[3]])
fields::image.plot(x, y, samples[[4]])
par(mfrow=c(1,1))




###############################################################################
# Exercise 8

source("include/GPFunctions2D.R")
source("include/PDEFunctions.R")
source("include/MonteCarloIntegration.R")
library(pracma) # used for interp2()
library(plotly) # Not required but makes nice plots you can drag around

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

# solves the PDE for many RFSamples. 
# Returns the solutions, mean and variance
RandomFieldMonteCarlo2D <- function(RFSamples, N=10, M=10, Lx=1, Ly=1, sigma2=1){
  # RFSamples:  list of Random Field samples (must include a f.base if required)
  # N:          x-direction summation limit
  # M:          y-direction summation limit
  # Lx:         x-axis domain, [0, Lx]
  # Ly:         y-axis domain, [0, Ly]
  
  no.samples <- length(RFSamples) # total number of samples
  RF.mean <- Reduce('+', RFSamples)/no.samples # mean of the samples
  
  # extract the dimensions of the problem from the samples
  I <- dim(RFSamples[[1]])[1]
  J <- dim(RFSamples[[1]])[2]
  
  # x and y-axis discretisations required for interpolation
  x <- seq(0, Lx, length=I)
  y <- seq(0, Ly, length=J)
  
  # build the empty data structures
  p.f <- lapply(1:no.samples, function(i) matrix(0, nrow=I, ncol=J))
  p.f.mean <- matrix(0, nrow=I, ncol=J)
  p.f.variance <- matrix(0, nrow=I, ncol=J)
  
  # iterate over the double sum
  for(n in 1:N){
    for(m in 1:M){
      
      #######
      # MEAN
      #######
      
      # estimate the spectral coefficients using Monte Carlo
      f.n.m <- MonteCarlo2D(user.func=function(x, y) {
        # interp2() linearly interpolates the matrix
        pracma::interp2(seq(0, Lx, length=I), seq(0, Ly, length=J), t(RF.mean), x, y, method="linear")*
          v.n.m(n, m, x, y)
      })$mean
      
      # use the coefficient to approximate the solution
      beta.n.m.evaluated <- beta.n.m(n, m, Lx=Lx, Ly=Ly)
      for(i in 1:length(x)){
        p.f.mean[i, ] <- p.f.mean[i, ] + f.n.m*v.n.m(n, m, x[i], y)/beta.n.m.evaluated
      }
      

      #######
      # Random Field samples
      #######
      
      # iterate over each samples
      for(k in 1:no.samples){
        # estimate the spectral coefficients using Monte Carlo
        f.n.m <- MonteCarlo2D(user.func=function(x, y) {
          interp2(seq(0, Lx, length=I), seq(0, Ly, length=J), t(RFSamples[[k]]), x, y, method="linear")*
            v.n.m(n, m, x, y)
        }, N=2000)$mean
        
        # use the coefficient to approximate the solution
        for(i in 1:length(x)){
          p.f[[k]][i, ] <- p.f[[k]][i, ] + f.n.m*v.n.m(n, m, x[i], y, Lx=Lx, Ly=Ly)/beta.n.m.evaluated
        }
      }
    }
  }
  
  # estimate the Monte Carlo variance from the p.f samples
  p.f.variance <- round(apply(array(unlist(p.f), c(I, J, N)), c(1, 2), var), 6)
  
  return(list("p.f"=p.f,
              "mean"=p.f.mean,
              "variance"=p.f.variance)
         )
}



I <- 30 # Number of points in x-direction
J <- I  # Number of points in y-direction
N <- 10 # number of samples

# This is the mean of random field (Gaussian process). Mentioned at 22:20 in Video2.mp4
f.base <- matrix(0.5, nrow=I, ncol=J)

# Build the x and y-axis discretisations
x <- seq(0, 1, length=I)
y <- seq(0, 1, length=J)

# Generate the Karhunen-Loeve Decomposition
# ONLY NEED TO RUN THIS ONCE
KLDecomp1 <- GetKLDecomposition(I=I, J=J, nu=1, sigma2=1, tau=1) # THIS IS VERY SLOW

# This is uses the decomposition above to generate new random samples and it very quick!
GP.samples <- GenerateSamples(N, KLDecomp1, I=I, J=J, f.base=f.base)

# plot the first random sample!
#fields::image.plot(x, y, GP.samples[[1]])

# solve the PDE
z <- GPMonteCarlo2D(GP.samples)

############
# solution

# Plot the mean solution
fig <- plot_ly(x=x, y=y, z=z$mean) %>% add_surface()
fig
#fields::image.plot(x, y, z$mean) # alternative plotting function

# Plot the first sample solution
fig <- plot_ly(x=x, y=y, z=z$p.f[[1]]) %>% add_surface()
fig
#fields::image.plot(x, y, z$p.f[[1]]) # alternative plotting function

############
# variance

fig <- plot_ly(x=x, y=y, z=z$variance) %>% add_surface()
fig
#fields::image.plot(x, y, z$variance) # alternative plotting function

#############
# histogram - of points (interpolated) at (0.5, 0.5)

central.values <- numeric(N)
for(n in 1:N){
  central.values[n] <- interp2(x, y, z$p.f[[n]], 0.5, 0.5, method="linear")
}
hist(central.values)


