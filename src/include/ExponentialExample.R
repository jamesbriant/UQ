u <- function(t, lambda, u0=2){
  # analytical system solution
  return(u0*exp(-t*lambda))
}

m <- function(t, lambda, epsilon, u0=2){
  # analytical mean
  a <- -u0/(2*epsilon*t)
  b <- exp(-t*(lambda+epsilon)) - exp(-t*(lambda-epsilon))
  
  return(a*b)
}