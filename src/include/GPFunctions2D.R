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
  
  return(A + t(A) + sigma2*diag(M))
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

GetKLDecomposition <- function(k=50, I=FALSE, J=FALSE, sigma2=1, nu=0.5, tau=1, Lx=1, Ly=1){
  # k:      Summation upper limit
  # I:      Number of cells along x-axis
  # J:      Number of cells along y-axis
  # sigma2: variance
  # nu:     smoothness variable
  # tau:    range variable
  # Lx:     x-axis domain, [0, Lx]
  # Ly:     y-axis domain, [0, Ly]
  
  ###########
  # CHECK INPUTS
  ###########
  
  if(I == FALSE && J == FALSE){
    print("ERROR: At least one of I or J must be set.")
    return(1)
  }else if(J == FALSE){
    J <- I
  }else if(I == FALSE){
    I <- J
  }
  
  M <- I*J
  
  ###########
  # IMPLEMENTATION
  ###########
  
  C <- GenerateC(M, I, sigma2=sigma2, nu=nu, tau=tau, Lx=Lx, Ly=Ly)
  eigen.sol <- SolveEigenProblem(C, k)
  
  return(eigen.sol)
}

GenerateSamples <- function(N, KLDecomp, f.base=FALSE, I=FALSE, J=FALSE){
  # N:        Number of 2D GP samples to draw
  # KLDecomp: returned decomposition from GetKLDecomposition()
  # f.base:   The base function (v-bar in the notes), must be pre-calculated and in matrix form of dimensions (I x J)
  # I:        Number of cells along x-axis
  # J:        Number of cells along y-axis
  
  if(I==FALSE && J==FALSE){
    I <- sqrt(length(KLDecomp$funcs[, 1]))
    J <- I
  }else if(I==FALSE){
    I <- length(KLDecomp$funcs[, 1])/J
  }else if(J==FALSE){
    J <- length(KLDecomp$funcs[, 1])/I
  }else if(length(KLDecomp$funcs[, 1]) != I*J){
    print("The dimensions requested do not agree with the decomposition provided.")
    return()
  }
  
  if(identical(f.base, FALSE)){
    f.base <- matrix(0, nrow=I, ncol=J)
  }else if(dim(f.base)[1] != I || dim(f.base)[2] != J){
    print("The dimensions of f.base do not agree with I and J.")
    return()
  }
  
  GPsamples <- list()
  
  for(n in 1:N){
    # initiate the nth sample
    GPsamples[[n]] <- f.base
    
    for(i in 1:KLDecomp$n){
      a <- sqrt(KLDecomp$values[i])*KLDecomp$funcs[, i]*rnorm(1, 0, 1)
      GPsamples[[n]] <- GPsamples[[n]] + matrix(a, byrow=FALSE, nrow=I, ncol=J)
    }
  }
  
  return(GPsamples)
}

GenerateGP2D <- function(N=1, k=50, I=FALSE, J=FALSE, sigma2=1, nu=0.5, tau=1, Lx=1, Ly=1, suppress.warning=FALSE){
  # N:      number of 2D GP samples to draw
  # k:      Summation upper limit
  # I:      Number of cells along x-axis
  # J:      Number of cells along y-axis
  # sigma2: variance
  # nu:     smoothness variable
  # tau:    range variable
  # Lx:     x-axis domain, [0, Lx]
  # Ly:     y-axis domain, [0, Ly]
  
  if(suppress.warning == FALSE){
    print("Avoid using this function. Instead, use GetKLDecomposition() and save the output. Then use GenerateSamples().")
  }
  
  # get the Karhunen-Loeve decomposition
  eigen.sol <- GetKLDecomposition(k, I, J, sigma2, nu, tau, Lx=1, Ly=1)
  
  return(GenerateSamples(N, eigen.sol, I=I, J=J))
}

FunctionToMatrix <- function(f, I, J, Lx=1, Ly=1){
  # f:      The base function (f in the notes), must accept 2 inputs (x,y)
  # I:      Number of cells along x-axis
  # J:      Number of cells along y-axis
  # Lx:     x-axis domain, [0, Lx]
  # Ly:     y-axis domain, [0, Ly]
  
  x <- seq(0, Lx, length=I)
  y <- seq(0, Lx, length=J)
  
  f.matrix <- matrix(0, nrow=I, ncol=J)
  
  for(j in 1:J){
    f.matrix[, j] <- sapply(x, f, y[j])
  }
  
  return(f.matrix)
}











