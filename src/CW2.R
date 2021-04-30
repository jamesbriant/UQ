library(pracma) # used for interp2()
source("include/MonteCarloIntegration.R")

GetWnm <- function(N, M=FALSE){
  # REQUIRES `include/MonteCarloIntegration.R` to be sourced!
  # REQUIRES `pracma` to be installed!

  raw.data <- read.delim("data/Data_groupB.txt", header=FALSE)
  x.grid <- unique(raw.data$V1)
  x.grid <- c(0, x.grid, 1)
  y.grid <- unique(raw.data$V2)
  y.grid <- c(0, y.grid, 1)
  p.eta <- matrix(0, nrow=length(x.grid), ncol=length(y.grid))
  for(i in 2:(length(x.grid)-1)){
    for(j in 2:(length(y.grid)-1)){
      p.eta[i, j] <- raw.data$V3[30*((i-1)-1) + (j-1)]
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
  
  #################
  # Find Wnm
  
  if(identical(M, FALSE)){
    M <- N
  }
  
  W.n.m <- matrix(0, nrow=N, ncol=M)
  
  for(n in 1:N){
    for(m in 1:M){
      # calculate the spectral coefficients
      W.n.m[n, m] <- cos(sqrt(beta.n.m(n, m))*0.322)*
        MonteCarlo2D(
          user.func=function(x, y){
            pracma::interp2(x.grid, y.grid, t(p.eta), x, y)*
            v.n.m(n, m, x, y)
          }, 
          N=5000
        )$mean
    }
  }
  
  return(W.n.m)
}


Test <- GetWnm(10, 10) # may take a little while to load...
fields::image.plot(1:10, 1:10, Test)













