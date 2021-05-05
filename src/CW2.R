library(pracma) # used for interp2()
#library(plotly)
source("include/MonteCarloIntegration.R")
source("include/GPFunctions2D.R")

raw.data <- read.delim("data/Data_groupB.txt", header=FALSE)
x.grid <- unique(raw.data$V1)
x.grid <- c(0, x.grid, 1)
y.grid <- unique(raw.data$V2)
y.grid <- c(0, y.grid, 1)
d.eta <- matrix(0, nrow=length(x.grid), ncol=length(y.grid))
for(i in 2:(length(x.grid)-1)){
  for(j in 2:(length(y.grid)-1)){
    d.eta[i, j] <- raw.data$V3[30*((i-1)-1) + (j-1)]
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
            pracma::interp2(x.grid, y.grid, t(d.eta), x, y)*
            v.n.m(n, m, x, y)
          }, 
          N=MonteCarloSize
        )$mean
    }
  }
  
  return(W.n.m)
}



N.max <- 18
M.max <- N.max
W.n.m <- GetWnm(N.max, x.grid, y.grid, d.eta, MonteCarloSize=20000) # may take a little while to load...


##########################################################
# Reconstruct the data using the Wnm decomposition

reconstruction <- matrix(0, nrow=length(x.grid), ncol=length(y.grid))

for(n in 1:N.max){
  for(m in 1:M.max){
    # evaluate the plotting points
    a <- 0.5*0.322*sqrt(beta.n.m(n, m, Lx=1, Ly=1))
    for(i in 1:length(x.grid)){
      reconstruction[i, ] <- reconstruction[i, ] + W.n.m[n, m]*cos(a)*v.n.m(n, m, x.grid[i], y.grid)
    }
  }
}

fields::image.plot(x.grid, y.grid, reconstruction, main=paste0("Reconstruction at t=t.star. N=", N.max))
fields::image.plot(x.grid, y.grid, d.eta, main="True Data")

#fig <- plot_ly(x=x.grid, y=y.grid, z=reconstruction) %>% add_surface()
#fig


###############################################################
# Find the initial condition - No UQ is done using this method
# Do not need to run this section

# 
# initial.condition <- matrix(0, nrow=length(x.grid), ncol=length(y.grid))
# 
# for(n in 1:N.max){
#   for(m in 1:M.max){
#     # evaluate the plotting points
#     for(i in 1:length(x.grid)){
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

# Set up the location of the centre of the cells
N.cells <- 50
F.grid.x <- seq(0, 1, length=N.cells+1)[1:N.cells]
h.x <- F.grid.x[2] - F.grid.x[1]
F.grid.x <- h.x/2 + F.grid.x
F.grid.y <- F.grid.x # y discretisation is the same as x discretisation
h.y <- h.x

F.h <- matrix(0, nrow=N.cells, ncol=N.cells) # matrix of Green's functions evaluations
data.interpolated <- matrix(0, nrow=N.cells, ncol=N.cells)
for(i in 1:N.cells){
  for(j in 1:N.cells){
    F.h[i, j] <- h.x*h.y*GreensFunction(F.grid.x[i], F.grid.y[i], 0, F.grid.x[j], F.grid.y[j])
    data.interpolated[i, j] <- pracma::interp2(x.grid, y.grid, t(d.eta), F.grid.x[i], F.grid.y[j])
  }
}

#fields::image.plot(F.grid.x, F.grid.y, F.h)

delta <- 1e-8

operator.inverted <- t(F.h) %*% inv(t(F.h) %*% F.h + delta*diag(N.cells))
solution <- operator.inverted %*% data.interpolated

fields::image.plot(F.grid.x, F.grid.y, solution)

#fig <- plot_ly(x=F.grid.x, y=F.grid.y, z=solution) %>% add_surface()
#fig




################################################
# This all needs flattening... gulp

f.prior <- matrix(0.5, nrow=N.cells, ncol=N.cells)
C.prior <- GenerateC(50^2, 50, sigma2=0.5, nu=0.5, tau=1) # THIS IS SLOW

A.temp <- C.prior %*% t(F.h) %*% inv(F.h %*% C.prior %*% t(F.h) + 0.01^2*diag(N.cells))
f.pos <- f.prior + A.temp %*% (data.interpolated - F.h %*% f.prior)
C.pos <- C.prior - A.temp %*% F.h %*% C.prior


fields::image.plot(F.grid.x, F.grid.y, f.pos)

fig <- plot_ly(x=F.grid.x, y=F.grid.y, z=f.pos) %>% add_surface()
fig














