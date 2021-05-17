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
  
  core.count <- detectCores()
  cl <- parallel::makeCluster(core.count-1)
  print(paste0("using ", core.count-1, " CPU cores."))
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
#F.h.0.322 <- GetF.h.parallel(0.322)
#F.h.0 <- GetF.h(0)
# F.h.0 <- GetF.h.parallel(0)
#F.h.0.4 <- GetF.h(0.4)
# F.h.0.4 <- GetF.h.parallel(0.4)

p.eta.vector <- matrix(0, nrow=x.grid.size*y.grid.size, ncol=1)
for(i in 1:x.grid.size){
  p.eta.vector[((i-1)*x.grid.size+1):(i*x.grid.size), 1] <- p.eta[i, ]
}

# This part below isn't needed, see Bayesian approach instead

# delta <- 1e-4
# 
# operator.inverted <- t(F.h.0.322) %*% inv(F.h.0.322 %*% t(F.h.0.322) + delta*diag(x.grid.size^2)) # equation 75 of notes
# solution.vector <- operator.inverted %*% p.eta.vector
# # now convert back to a matrix
# solution <- matrix(solution.vector, nrow=F.x.grid.size, ncol=F.y.grid.size, byrow=TRUE)
# 
# fields::image.plot(F.x.grid, F.y.grid, solution)

#fig <- plot_ly(x=F.x.grid, y=F.y.grid, z=solution) %>% add_surface()
#fig




################################################
# These are the Bayesian steps

# prior is GP with flat mean and Matern covariance
f.prior <- matrix(0.2, nrow=N.cells^2, ncol=1)
C.prior <- GenerateC(N.cells^2, N.cells, sigma2=0.1, nu=2, tau=0.1) # THIS IS SLOW, ONLY RUN ME ONCE

A.temp <- C.prior %*% t(F.h.0.322) %*% inv(F.h.0.322 %*% C.prior %*% t(F.h.0.322) + 0.01^2*diag(x.grid.size*y.grid.size))
f.pos <- f.prior + A.temp %*% (p.eta.vector - F.h.0.322 %*% f.prior)
f.pos.matrix <- matrix(f.pos, nrow=N.cells, byrow=TRUE)
C.pos <- C.prior - A.temp %*% F.h.0.322 %*% C.prior
variance.matrix <- matrix(diag(C.pos), nrow=N.cells, byrow=TRUE)

# this is the posterior distribution
fields::image.plot(F.x.grid, F.y.grid, f.pos.matrix)

#fig <- plot_ly(x=F.x.grid, y=F.y.grid, z=f.pos.matrix) %>% add_surface()
#fig


# This is the variance of the posterior distribution
fields::image.plot(F.x.grid, F.y.grid, sqrt(variance.matrix))

#fig <- plot_ly(x=F.x.grid, y=F.y.grid, z=sqrt(variance.matrix)) %>% add_surface()
#fig




################################################################################
# Finding the displacement at Lego Town

#######################
# Method 1
# This is not used due to the random nature of the estimates of the Wnm terms.

# GetWnm <- function(N, x.grid, y.grid, data.matrix, c=0.5, M=FALSE, MonteCarloSize=5000, silent=TRUE){
#   # REQUIRES `include/MonteCarloIntegration.R` to be sourced!
#   # REQUIRES `pracma` to be installed!
#   
#   #################
#   # Find Wnm
#   
#   if(identical(M, FALSE)){
#     M = N
#   }
#   
#   W.n.m <- matrix(0, nrow=N, ncol=M)
#   
#   for(n in 1:N){
#     if(silent==FALSE){
#       print(n)
#     }
#     for(m in 1:M){
#       # calculate the spectral coefficients
#       W.n.m[n, m] <- MonteCarlo2D(
#         user.func=function(x, y){
#           pracma::interp2(x.grid, y.grid, t(data.matrix), x, y)*
#             v.n.m(n, m, x, y)
#         }, 
#         N=MonteCarloSize
#       )$mean
#     }
#   }
#   
#   return(W.n.m)
# }
# 
# GetTimeToLegoCity <- function(x.start, y.start, c=0.5){
#   d <- sqrt((0.495 - x.start)^2 + (0.495 - y.start)^2)
#   return(d/c)
# }
# 
# MasterBuilderEmmet <- function(W, initial.condition, F.x.grid, F.y.grid, c=0.5){ # I hope you've seen the Lego Movie, otherwise this function name won't mean anything :D
#   # This function builds the system state whenever the wave is expected to hit Lego City
#   
#   # Find the epicentre of the earthquake
#   centre.location.encoded <- which(initial.condition==max(initial.condition))
#   centre.location.decoded <- DecodeLocation2D(centre.location.encoded, length(F.x.grid))
#   centre.location.x <- F.x.grid[centre.location.decoded[1]]
#   centre.location.y <- F.y.grid[centre.location.decoded[2]]
#   
#   # Find the time the epicentre reaches Lego City
#   t <- GetTimeToLegoCity(centre.location.x, centre.location.y)
#   
#   # build the entire domain
#   output <- matrix(0, nrow=length(F.x.grid), ncol=length(F.y.grid))
#   for(n in 1:dim(W)[1]){
#     for(m in 1:dim(W)[2]){
#       for(i in 1:length(F.x.grid)){
#         output[i, ] <- output[i, ] + W[n, m]*v.n.m(n, m, F.x.grid[i], F.y.grid)*cos(pi*c*t*sqrt(n^2 + m^2))
#       }
#     }
#   }
#   
#   # return the displacement in the domain
#   return(output)
# }
# 
# MasterBuilderEmmet2 <- function(W, initial.condition, F.x.grid, F.y.grid, c=0.5){ # I hope you've seen the Lego Movie, otherwise this function name won't mean anything :D
#   # This is an improved version of MasterBuilderEmmet(), but only returns the value at Lego City.
#   # The rest of the domain is not calculated
#   
#   # Find the epicentre of the earthquake
#   centre.location.encoded <- which(initial.condition==max(initial.condition))
#   centre.location.decoded <- DecodeLocation2D(centre.location.encoded, length(F.x.grid))
#   centre.location.x <- F.x.grid[centre.location.decoded[1]]
#   centre.location.y <- F.y.grid[centre.location.decoded[2]]
#   
#   # Find the time the epicentre reaches Lego City
#   t <- GetTimeToLegoCity(centre.location.x, centre.location.y)
#   
#   # Find the F.x.grid locations surrounding Lego City. We interpolate later
#   x.a <- max(which(F.x.grid < 0.495))
#   x.b <- min(which(F.x.grid > 0.495))
#   y.a <- max(which(F.y.grid < 0.495))
#   y.b <- min(which(F.y.grid > 0.495))
#   
#   # build the region of interest
#   output <- matrix(0, nrow=2, ncol=2)
#   for(n in 1:dim(W)[1]){
#     for(m in 1:dim(W)[2]){
#       temp <- W[n, m]*cos(pi*c*t*sqrt(n^2 + m^2))
#       output[1, ] <- output[1, ] + temp*v.n.m(n, m, F.x.grid[x.a], F.y.grid[y.a:y.b])
#       output[2, ] <- output[2, ] + temp*v.n.m(n, m, F.x.grid[x.b], F.y.grid[y.a:y.b])
#     }
#   }
#   
#   # interpolate the displacement at Lego City and return
#   return(pracma::interp2(F.x.grid[x.a:x.b], F.y.grid[y.a:y.b], t(output), 0.495, 0.495))
# }
# 
# # add on the border for interpolation later
# F.x.grid.border <- c(0, F.x.grid, 1)
# F.y.grid.border <- c(0, F.y.grid, 1)
# f.pos.matrix.0.border <- matrix(0, nrow=F.x.grid.size+2, ncol=F.y.grid.size+2)
# f.pos.matrix.0.border[2:(F.x.grid.size+1), 2:(F.x.grid.size+1)] <- f.pos.matrix
# 
# N.max <- 18 # check this value!
# M.max <- N.max
# W.n.m <- GetWnm(N.max, 
#                 F.x.grid.border, 
#                 F.y.grid.border, 
#                 f.pos.matrix.0.border, 
#                 MonteCarloSize=10000,
#                 silent=FALSE) # may take a little while to load...
# 
# 
# output <- MasterBuilderEmmet(W.n.m, f.pos.matrix.0.border, F.x.grid.border, F.y.grid.border)
# fields::image.plot(F.x.grid.border, F.y.grid.border, output)
# #pracma::interp2(F.x.grid.border, F.y.grid.border, t(output), 0.495, 0.495) # used to check the output is the same as MasterBuilderEmmet2(). It was the same!
# 
# MasterBuilderEmmet2(W.n.m, f.pos.matrix.0.border, F.x.grid.border, F.y.grid.border)



########################
# Method 2
# Propagate forwards samples from N(f.pos, C.pos) using F.h.forward

GetTimeToLegoCity <- function(x.start, y.start, c=0.5){
  d <- sqrt((0.495 - x.start)^2 + (0.495 - y.start)^2)
  return(d/c)
}

GetF.h.forward <- function(t){
  F.h <- matrix(0, nrow=(N.cells)^2, ncol=(N.cells)^2) # matrix of Green's functions evaluations
  pb <- txtProgressBar(min = 0, max = (N.cells)^2, style = 3) # progress bar
  for(i in 1:((N.cells)^2)){
    # clever matrix location decoder
    temp <- DecodeLocation2D(i, (N.cells))
    F.x.grid.i <- temp[1]
    F.y.grid.i <- temp[2]
    
    for(j in 1:(N.cells)){
      F.h[i, ((j-1)*(N.cells)+1):(j*(N.cells))] <- h.x*h.y*GreensFunction(F.x.grid[F.x.grid.i], F.y.grid[F.y.grid.i], t, F.x.grid[j], F.y.grid)
    }
    
    setTxtProgressBar(pb, i) # update progress bar
  }
  
  return(F.h)
}

GetF.h.forward.parallel <- function(t){
  library(foreach)
  library(parallel)
  library(doParallel)
  
  core.count <- detectCores()
  cl <- parallel::makeCluster(core.count-1)
  print(paste0("using ", core.count-1, " CPU cores."))
  doParallel::registerDoParallel(cl)
  F.h <- foreach(i = 1:((N.cells)^2), .combine='rbind', .export = ls(globalenv())) %dopar% {
    # clever matrix location decoder
    temp <- DecodeLocation2D(i, (N.cells))
    F.x.grid.i <- temp[1]
    F.y.grid.i <- temp[2]
    
    output <- numeric((N.cells)^2)
    
    for(j in 1:(N.cells)){
      output[((j-1)*(N.cells)+1):(j*(N.cells))] <- h.x*h.y*GreensFunction(F.x.grid[F.x.grid.i], F.y.grid[F.y.grid.i], t, F.x.grid[j], F.y.grid)
    }
    
    output
  }
  parallel::stopCluster(cl)
  
  return(F.h)
}

# only run these 3 lines once.
KL.decomp <- SolveEigenProblem(C.pos) # this is very slow
saved.F.hs.times <- numeric()
saved.F.hs <- list()

CalculateF.hs <- function(pos.samples, saved.F.hs, saved.F.hs.times, parallel=FALSE){
  for(i in 1:length(pos.samples)){
    # Find the epicentre of the earthquake
    centre.location.encoded <- which(pos.samples[[i]]==max(pos.samples[[i]]))
    centre.location.decoded <- DecodeLocation2D(centre.location.encoded, N.cells)
    centre.location.x <- F.x.grid[centre.location.decoded[1]]
    centre.location.y <- F.y.grid[centre.location.decoded[2]]
    
    # Find the time the epicentre reaches Lego City
    t <- GetTimeToLegoCity(centre.location.x, centre.location.y)
    
    if(!(t %in% saved.F.hs.times)){
      print(paste0("New time ", t))
      
      saved.F.hs.times <- c(saved.F.hs.times, t)
      
      if(parallel==FALSE){
        saved.F.hs[[length(saved.F.hs.times)]] <- GetF.h.forward(t)
      }else{
        saved.F.hs[[length(saved.F.hs.times)]] <- GetF.h.forward.parallel(t)
      }
    }
  }
  
  return(list("saved.F.hs"=saved.F.hs,
              "saved.F.hs.times"=saved.F.hs.times))
}




sample.count <- 10000
pos.samples <- GenerateSamples(sample.count, KL.decomp, f.base=f.pos.matrix)

run.F.h.search <- CalculateF.hs(pos.samples, saved.F.hs, saved.F.hs.times, parallel=FALSE)
saved.F.hs.times <- run.F.h.search$saved.F.hs.times
saved.F.hs <- run.F.h.search$saved.F.hs

lego.city.displacements <- numeric(sample.count)
pb <- txtProgressBar(min = 0, max = sample.count, style = 3) # progress bar
for(i in 1:sample.count){
  # Find the epicentre of the earthquake
  centre.location.encoded <- which(pos.samples[[i]]==max(pos.samples[[i]]))
  centre.location.decoded <- DecodeLocation2D(centre.location.encoded, N.cells)
  centre.location.x <- F.x.grid[centre.location.decoded[1]]
  centre.location.y <- F.y.grid[centre.location.decoded[2]]
  
  # Find the time the epicentre reaches Lego City
  t <- GetTimeToLegoCity(centre.location.x, centre.location.y)
  
  # find the location of the right F.h(t)
  F.h.time <- which(saved.F.hs.times == t)
  
  # DO NOT USE F.h.forward HERE. THIS FUNCTION SHOULD BE SPECIFIC TO EACH SAMPLE
  u.lego.city.i <- matrix(saved.F.hs[[F.h.time]] %*% matrix(pos.samples[[i]], ncol=1), nrow=N.cells, byrow=TRUE)
  
  # interpolate the location for lego city
  lego.city.displacements[i] <- pracma::interp2(F.x.grid, F.y.grid, t(u.lego.city.i), 0.495, 0.495)

  setTxtProgressBar(pb, i) # update progress bar
}


hist(lego.city.displacements, freq=FALSE, main="Histogram of Maximum Lego Town Displacements from 10000 Samples", xlab="Maximum Displacement")
lines(seq(0, 0.2, length=200), dnorm(seq(0, 0.2, length=200), mean(lego.city.displacements), sd=sqrt(var(lego.city.displacements))), col="blue")
abline(v=0.12, col="red")

# probability of buildings collapsing
sum(lego.city.displacements >= 0.12)/sample.count 




######################################
# Method 3
# This doesn't work, I don't know why. Stick with method 2

# F.h.0.422 <- GetF.h(0.422)
# u.0.1.method3.vector <- t(F.h.0.422) %*% p.eta.vector
# u.0.1.method3 <- matrix(u.0.1.method3.vector, nrow=N.cells, byrow=TRUE)
# fields::image.plot(F.x.grid, F.y.grid, u.0.1.method3)


