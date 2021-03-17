#############################################################################
# SET THIS FOLDER AS THE WORKING DIRECTORY
#############################################################################

source("include/Quadrature.R")
source("include/DiscreteProjections.R")
source("include/ClassTestFunctions.R")
source("include/QMC.R")

#############################################################################
# QUESTION a)

# Use the Wigner-Chebyshev scheme 

lambda <- seq(-1, 1, length=300)
plot(lambda, u(lambda), type="l", ylim=c(-0.2, 1.5), main="Wigner-Chebyshev Scheme of different degree")

Js <- c(5, 10, 20, 30)

coefficients.u <- list()
gammas.u <- list()
u.n <- list()

# This loop takes a while on my laptop...
for(i in 1:length(Js)){
  sol <- WignerChebyshev(lambda, 20, u, Js[i])
  
  # plot the output
  lines(lambda, sol$approximation, col=i+1)
  
  # save the coefficients
  coefficients.u[[i]] <-  sol$coefficients
  gammas.u[[i]] <- sol$gammas
}

legend("topright", 
       legend=c("u(lambda)", sapply(1:length(Js), function(i) paste0("J=", Js[i]))),
       col=1:(length(Js)+1),
       lty=1
       )

# I think something is wrong with the scale of my projections? I am not sure what. :(

#############################################################################
# QUESTION b)

means.u <- numeric(length(Js))
vars.u <- numeric(length(Js))
for(i in 1:length(Js)){
  means.u[i] <- coefficients.u[[i]][1]
  vars.u[i] <- sum(coefficients.u[[i]]^2*gammas.u[[i]])
}

# The Monte Carlo mean and variance is
n <- 10^7
wigner.sample <- rwigner(n)
u.MC <- u(wigner.sample)
mean.MC <- mean(u.MC)

var.MC <- mean(u.MC^2) - mean.MC^2



#The means and variances for J=5, 10, 20 and 30 are
means.u
vars.u

# MC mean and variance
mean.MC
var.MC

# The means for the MC estimates and the spectral estimates are similar. 
# But the variance terms are very different. This could be partly blamed on Runge's phenomenon.

#############################################################################
# QUESTION c)

wigner.sample.2 <- rwigner(10^4)
hist(u(wigner.sample.2))


sample.5 <- numeric(10^4)
for(i in 0:5){
  sample.5 <- sample.5 + coefficients.u[[1]][i] * ChebyshevPolynomial(i, wigner.sample.2)
}


hist(sample.5)

sample.20 <- numeric(10^4)
for(i in 0:20){
  sample.20 <- sample.20 + coefficients.u[[4]][i] * ChebyshevPolynomial(i, wigner.sample.2)
}


hist(sample.20)
