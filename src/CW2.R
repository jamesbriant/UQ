library(pracma) # used for interp2()
#library(plotly)
source("include/MonteCarloIntegration.R")

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



N.max <- 18
M.max <- N.max
W.n.m <- GetWnm(N.max, x.grid, y.grid, p.eta, N=20000) # may take a little while to load...


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
fields::image.plot(x.grid, y.grid, p.eta, main="True Data")

#fig <- plot_ly(x=x.grid, y=y.grid, z=reconstruction) %>% add_surface()
#fig


###############################################################
# Find the initial condition


initial.condition <- matrix(0, nrow=length(x.grid), ncol=length(y.grid))

for(n in 1:N.max){
  for(m in 1:M.max){
    # evaluate the plotting points
    for(i in 1:length(x.grid)){
      initial.condition[i, ] <- initial.condition[i, ] + W.n.m[n, m]*v.n.m(n, m, x.grid[i], y.grid)
    }
  }
}

fields::image.plot(x.grid, y.grid, initial.condition, main=paste0("Initial Condition. N=", N.max))

#fig <- plot_ly(x=x.grid, y=y.grid, z=initial.condition) %>% add_surface()
#fig































