#############
# LOADING PACKAGES
#############

LoadRSpectra <- function(){
  if(!("RSpectra" %in% (.packages()))){
    library(RSpectra)
  }
}

#############
# LOCATION ENCODING & DECODING
#############

EncodeLocation2D <- function(i, j, I){
  return((i-1)*I + j)
}

DecodeLocation2D <- function(m, I){
  i <- ((m-1)%%I) + 1
  j <- ceiling((m-0.5)/I)
  
  return(c(i, j))
}

# PhiBasis2D <- function(m, eta, I, Lx=1, Ly=1){
#   # m:    encoded basis number
#   # eta:  2D location vector
#   # Lx:   x-axis domain, [0, Lx]
#   # Ly:   y-axis domain, [0, Ly]
#   
#   m.decoded <- DecodeLocation2D(m, I)
#   i <- m.decoded[1]
#   j <- m.decoded[2]
#   
#   a <- 2/sqrt(Lx*Ly)
#   b <- sin(i*pi*eta[1]/Lx)
#   c <- sin(j*pi*eta[2]/Ly)
#   
#   return(a*b*c)
# }

#############
# GENERAL FUNCTIONS
#############

Matern2D <- function(x1, x2, sigma2=1, nu=0.5, tau=1){
  # x1:     2D vector
  # x2:     2D vector
  # sigma2: variance
  # nu:     smoothness variable
  # tau:    range variable
  
  d <- sqrt((x1[1] - x2[1])^2 + (x1[2] - x2[2])^2)/tau
  
  a <- sigma2/(2^(nu-1)*gamma(nu))
  b <- (d)^nu
  c <- besselK(d, nu)
  
  return(a*b*c)
}

CalculateEta <- function(m, I, J, Lx=1, Ly=1){
  # i:  encoded location
  # I:  Number of cells along x-axis
  # J:  Number of cells along y-axis
  # Lx: x-axis domain, [0, Lx]
  # Ly: y-axis domain, [0, Ly]
  
  m.decoded <- DecodeLocation2D(m, I)
  i <- m.decoded[1]
  j <- m.decoded[2]
  
  h.x <- Lx/(I+1)
  h.y <- Ly/(J+1)
  
  x.dot <- h.x*(i - 0.5)
  y.dot <- h.y*(j - 0.5)
  
  return(c(x.dot, y.dot))
}

GenerateC <- function(M, I, sigma2=1, nu=0.5, tau=1, Lx=1, Ly=1){
  # M:      Total number of evaluated points, I*J
  # I:      Number of cells along x-axis
  # sigma2: variance
  # nu:     smoothness variable
  # tau:    range variable
  # Lx:     x-axis domain, [0, Lx]
  # Ly:     y-axis domain, [0, Ly]
  
  # Since C is symmetric, only calculate lower triangular and add transpose.
  # Leading diagonal is 0 everywhere (as d = ||x-x|| = 0).
  
  A <- matrix(0, ncol=M, nrow=M)
  
  # calculate the number of cells in y-direction
  J <- M/I
  
  # fill in the lower-triangle
  for(i in 2:M){
    for(j in 1:(i-1)){
      eta.i <- CalculateEta(i, I, J, Lx=1, Ly=1)
      eta.j <- CalculateEta(j, I, J, Lx=1, Ly=1)
      
      A[i, j] <- Matern2D(eta.i, eta.j, sigma2=sigma2, nu=nu, tau=tau)
    }
  }
  
  return(A+t(A))
}

SolveEigenProblem <- function(C, k=50){
  # C: a square matrix
  # k: the first k eigenvalues and eigenvectors will be returned
  
  LoadRSpectra()
  
  sol <- eigs_sym(C, k, which="LM")
  
  n <- max(which(sol$values > 0))
  values <- sol$values[1:n]
  funcs <- sol$vectors[, 1:n]
  
  return(list(n=n,
              values=values,
              funcs=funcs)
         )
}

GenerateGP2D <- function(N=1, k=50, M=FALSE, I=FALSE, sigma2=1, nu=0.5, tau=1, Lx=1, Ly=1){
  # N:      number of 2D GP samples to draw
  # k:      Summation upper limit
  # M:      Total number of evaluated points, I*J
  # I:      Number of cells along x-axis
  # sigma2: variance
  # nu:     smoothness variable
  # tau:    range variable
  # Lx:     x-axis domain, [0, Lx]
  # Ly:     y-axis domain, [0, Ly]
  
  ###########
  # CHECK INPUTS
  ###########
  
  if(M == FALSE && I == FALSE){
    print("ERROR: At least one of M or I must be set.")
    return(1)
  }else if(M == FALSE){
    M <- I^2
  }else if(I == FALSE){
    I <- sqrt(M)
  }else if(M %% I != 0){
    print("ERROR: I must divide M exactly.")
    return(2)
  }
  
  ###########
  # IMPLEMENTATION
  ###########
  
  C <- GenerateC(M, I)
  eigen.sol <- SolveEigenProblem(C, k)
  
  J <- M/I
  
  GPsamples <- list()
  
  for(n in 1:N){
    # initiate the nth sample
    GPsamples[[n]] <- matrix(0, nrow=I, ncol=J)
    
    for(i in 1:eigen.sol$n){
      a <- sqrt(eigen.sol$values[i])*eigen.sol$funcs[, i]*rnorm(1, 0, 1)
      GPsamples[[n]] <- GPsamples[[n]] + matrix(a, byrow=FALSE, nrow=I, ncol=J)
    }
    
    #GPsamples[[n]] <- t(GPsamples[[n]])
  }
  
  return(GPsamples)
}















