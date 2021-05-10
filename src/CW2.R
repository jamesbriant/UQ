library(pracma) # used for interp2()
#library(plotly)
source("include/MonteCarloIntegration.R")
source("include/GPFunctions2D.R")

raw.data <- read.delim("data/Data_groupB.txt", header=FALSE)
x.grid <- unique(raw.data$V1)
#x.grid <- c(0, x.grid, 1)
x.grid.size <- length(x.grid)
y.grid <- unique(raw.data$V2)
#y.grid <- c(0, y.grid, 1)
y.grid.size <- length(y.grid)
p.eta <- matrix(0, nrow=x.grid.size, ncol=y.grid.size)
if(x.grid.size == 32){
  for(i in 2:(x.grid.size-1)){
    for(j in 2:(y.grid.size-1)){
      p.eta[i, j] <- raw.data$V3[30*((i-1)-1) + (j-1)]
    }
  }
}else{
  for(i in 1:x.grid.size){
    for(j in 1:y.grid.size){
      p.eta[i, j] <- raw.data$V3[30*(i-1) + (j)]
    }
  }
}

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

GetWnm <- function(N, x.grid, y.grid, data.matrix, c=0.5, M=FALSE, MonteCarloSize=5000){
  # REQUIRES `include/MonteCarloIntegration.R` to be sourced!
  # REQUIRES `pracma` to be installed!

  #################
  # Find Wnm
  
  if(identical(M, FALSE)){
    M = N
  }

  W.n.m <- matrix(0, nrow=N, ncol=M)
  
  for(n in 1:N){
    for(m in 1:M){
      # calculate the spectral coefficients
      W.n.m[n, m] <- 1/cos(sqrt(beta.n.m(n, m))*0.322*c)*
        MonteCarlo2D(
          user.func=function(x, y){
            pracma::interp2(x.grid, y.grid, t(p.eta), x, y)*
            v.n.m(n, m, x, y)
          }, 
          N=MonteCarloSize
        )$mean
    }
  }
  
  return(W.n.m)
}



# N.max <- 18
# M.max <- N.max
# W.n.m <- GetWnm(N.max, x.grid, y.grid, p.eta, MonteCarloSize=20000) # may take a little while to load...


##########################################################
# Reconstruct the data using the Wnm decomposition

# reconstruction <- matrix(0, nrow=x.grid.size, ncol=y.grid.size)
# 
# for(n in 1:N.max){
#   for(m in 1:M.max){
#     # evaluate the plotting points
#     a <- 0.5*0.322*sqrt(beta.n.m(n, m, Lx=1, Ly=1))
#     for(i in 1:x.grid.size){
#       reconstruction[i, ] <- reconstruction[i, ] + W.n.m[n, m]*cos(a)*v.n.m(n, m, x.grid[i], y.grid)
#     }
#   }
# }
# 
# fields::image.plot(x.grid, y.grid, reconstruction, main=paste0("Reconstruction at t=t.star. N=", N.max))
# fields::image.plot(x.grid, y.grid, p.eta, main="True Data")

#fig <- plot_ly(x=x.grid, y=y.grid, z=reconstruction) %>% add_surface()
#fig


###############################################################
# Find the initial condition - No UQ is done using this method
# Do not need to run this section

# 
# initial.condition <- matrix(0, nrow=x.grid.size, ncol=y.grid.size)
# 
# for(n in 1:N.max){
#   for(m in 1:M.max){
#     # evaluate the plotting points
#     for(i in 1:x.grid.size){
#       initial.condition[i, ] <- initial.condition[i, ] + W.n.m[n, m]*v.n.m(n, m, x.grid[i], y.grid)
#     }
#   }
# }
# 
# fields::image.plot(x.grid, y.grid, initial.condition, main=paste0("Initial Condition. N=", N.max))
# 
# #fig <- plot_ly(x=x.grid, y=y.grid, z=initial.condition) %>% add_surface()
# #fig



###########################################################
#

GreensFunction <- function(x1, y1, t, x2, y2, N.max=18, M.max=FALSE, c=0.5){
  if(identical(M.max, FALSE)){
    M.max <- N.max
  }
  
  gf <- 0
  
  for(n in 1:N.max){
    for(m in 1:M.max){
      gf <- gf + v.n.m(n, m, x1, y1)*v.n.m(n, m, x2, y2)*cos(c*t*pi*sqrt(n^2 + m^2))
    }
  }
  
  return(gf)
}

##############
# Set up the location of the centre of the cells
N.cells <- 50 # Number of cells in x and y-direction

# Set up the x-direction discretisation
F.x.grid <- seq(0, 1, length=N.cells+1)[1:N.cells]
h.x <- F.x.grid[2] - F.x.grid[1]
F.x.grid <- h.x/2 + F.x.grid
F.x.grid.size <- length(F.x.grid)

# Set up the y-direction discretisation
F.y.grid <- F.x.grid # y discretisation is the same as x discretisation
F.y.grid.size <- length(F.y.grid)
h.y <- h.x


GetF.h <- function(t){
  F.h <- matrix(0, nrow=x.grid.size^2, ncol=N.cells^2) # matrix of Green's functions evaluations
  pb <- txtProgressBar(min = 0, max = x.grid.size^2, style = 3) # progress bar
  for(i in 1:(x.grid.size*y.grid.size)){
    # clever matrix location decoder
    temp <- DecodeLocation2D(i, x.grid.size)
    x.grid.i <- temp[1]
    y.grid.i <- temp[2]
    
    for(j in 1:N.cells){
      F.h[i, ((j-1)*N.cells+1):(j*N.cells)] <- h.x*h.y*GreensFunction(x.grid[x.grid.i], y.grid[y.grid.i], t, F.x.grid[j], F.y.grid)
    }
    
    setTxtProgressBar(pb, i) # update progress bar
  }
  
  return(F.h)
}

GetF.h.parallel <- function(t){
  library(foreach)
  library(parallel)
  library(doParallel)
  
  cl <- parallel::makeCluster(detectCores()-1)
  doParallel::registerDoParallel(cl)
  F.h <- foreach(i = 1:(x.grid.size*y.grid.size), .combine='rbind', .export = ls(globalenv())) %dopar% {
    # clever matrix location decoder
    temp <- DecodeLocation2D(i, x.grid.size)
    x.grid.i <- temp[1]
    y.grid.i <- temp[2]

    output <- numeric(N.cells^2)

    for(j in 1:N.cells){
      output[((j-1)*N.cells+1):(j*N.cells)] <- h.x*h.y*GreensFunction(x.grid[x.grid.i], y.grid[y.grid.i], t, F.x.grid[j], F.y.grid)
    }

    output
  }
  parallel::stopCluster(cl)
  
  return(F.h)
}

# Generate the F^h matrix in the notes
F.h.0.322 <- GetF.h(0.322)
# F.h.0.322 <- GetF.h.parallel(0.322)
#F.h.0 <- GetF.h(0)
# F.h.0 <- GetF.h.parallel(0)
#F.h.0.4 <- GetF.h(0.4)
# F.h.0.4 <- GetF.h.parallel(0.4)

p.eta.vector <- matrix(0, nrow=x.grid.size*y.grid.size, ncol=1)
for(i in 1:x.grid.size){
  p.eta.vector[((i-1)*x.grid.size+1):(i*x.grid.size), 1] <- p.eta[i, ]
}

delta <- 1e-4

operator.inverted <- t(F.h.0.322) %*% inv(F.h.0.322 %*% t(F.h.0.322) + delta*diag(x.grid.size^2)) # equation 75 of notes
solution.vector <- operator.inverted %*% p.eta.vector
# now convert back to a matrix
solution <- matrix(solution.vector, nrow=F.x.grid.size, ncol=F.y.grid.size, byrow=TRUE)

fields::image.plot(F.x.grid, F.y.grid, solution)

#fig <- plot_ly(x=F.x.grid, y=F.y.grid, z=solution) %>% add_surface()
#fig




################################################
# These are the Bayesian steps

# prior is GP with flat mean and Matern covariance
f.prior <- matrix(0.5, nrow=N.cells^2, ncol=1)
C.prior <- GenerateC(N.cells^2, N.cells, sigma2=0.1, nu=1.5, tau=1) # THIS IS SLOW, ONLY RUN ME ONCE

A.temp <- C.prior %*% t(F.h.0.322) %*% inv(F.h.0.322 %*% C.prior %*% t(F.h.0.322) + 0.01^2*diag(x.grid.size*y.grid.size))
f.pos <- f.prior + A.temp %*% (p.eta.vector - F.h.0.322 %*% f.prior)
f.pos.matrix <- matrix(f.pos, nrow=N.cells, byrow=TRUE)
C.pos <- C.prior - A.temp %*% F.h.0.322 %*% C.prior
variance.matrix <- matrix(diag(C.pos), nrow=N.cells, byrow=TRUE)

# this is the posterior distribution
fields::image.plot(F.x.grid, F.y.grid, f.pos.matrix)

fig <- plot_ly(x=F.x.grid, y=F.y.grid, z=f.pos.matrix) %>% add_surface()
fig


# This is the variance of the posterior distribution
# These values should be positive... what has gone wrong? :(
fields::image.plot(F.x.grid, F.y.grid, sqrt(-variance.matrix))

fig <- plot_ly(x=F.x.grid, y=F.y.grid, z=sqrt(-variance.matrix)) %>% add_surface()
fig












