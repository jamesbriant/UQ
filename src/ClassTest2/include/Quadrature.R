LoadOrthopolynom <- function(){
  if(!("orthopolynom" %in% (.packages()))){
    library(orthopolynom)
  }
}

CompositeTrapezoidalRule <- function(user.func, N=5, a=0, b=1){
  # N is the number of nodes
  
  h <- (b-a)/(N-1)
  x <- seq(a, b, by=h)
  
  quad.points <- sapply(x, user.func)
  
  approximation <- h*(quad.points[1] + 2*quad.points[2:(N-1)] + quad.points[N])/2
  
  return(approximation)
}

SimpsonRule <- function(user.func, a=0, b=1){
  mid.point <- (a+b)/2
  
  return((b-a)*(user.func(a) + 4*user.func(mid.point) + user.func(b))/6)
}

CompositeSimpsonRule <- function(user.func, N=5, a=0, b=1){
  # N is the number of nodes
  
  h <- (b-a)/(N-1)
  x <- seq(a, b, by=h)
  
  approximation <- 0
  
  for(i in 1:(length(x)-1)){
    approximation <- approximation + SimpsonRule(user.func, a=x[i], b=x[i+1])
  }
  
  return(approximation)
}

LegendrePolynomial <- function(n, x){
  if(n == 1) return(x)
  if(n == 0) return(1)
  
  return((2*n - 1)*x*LegendrePolynomial(n-1, x)/n - (n-1)*LegendrePolynomial(n-2, x)/n)
}

LegendrePolynomialDerivative <- function(n, x){
  if(n == 1) return(1)
  if(n == 0) return(0)
  
  a <- (2*n - 1)*LegendrePolynomial(n-1, x)/n
  b <- (2*n - 1)*x*LegendrePolynomialDerivative(n-1, x)/n
  c <- (n - 1)*LegendrePolynomialDerivative(n-2, x)/n
  
  return(a + b + c)
}

GaussQuadrature <- function(user.func, J, lambda.dagger=0, epsilon=1){
  # J is the polynomial degree
  # lambda.dagger and epsilon map [-1, 1] on to the region of interest
  
  LoadOrthopolynom()
  
  r <- legendre.recurrences(J, normalized=TRUE)
  m.r <- monic.polynomial.recurrences(r)
  roots <- polynomial.roots(m.r)[[J+1]]
  
  lambda <- lambda.dagger + epsilon*roots
  
  #approximation <- 0
  
  #for(i in 1:J){
  #  weight <- 2/(LegendrePolynomial(J-1, roots[i])*LegendrePolynomialDerivative(J, roots[i])*J)
  #  approximation <- approximation + weight*user.func(lambda.dagger + epsilon*roots[i])
  #}
  
  weights <- 2/(LegendrePolynomial(J-1, roots)*LegendrePolynomialDerivative(J, roots)*J)
  approximation <- sum(weights * user.func(lambda))
  
  return(approximation)
}

HermitePolynomial <- function(n, x){
  if(n == 1) return(x)
  if(n == 0) return(1)
  
  return(x*HermitePolynomial(n-1, x) - (n-1)*HermitePolynomial(n-2, x))
}

LaguerrePolynomial <- function(n, x){
  if(n == 1) return(1-x)
  if(n == 0) return(1)
  
  return(((2*n - 1 - x)*LaguerrePolynomial(n-1, x) - (n-1)*LaguerrePolynomial(n-2, x))/n)
}

JacobiPolynomial <- function(n, a, b, x){
  if(n == 1) return((a+1) + (a+b+2)*(x-1)/2)
  if(n == 0) return(1)
  
  denom <- 2*n*(n+a+b)*(2*n+a+b-2)
  part1 <- (2*n+a+b-1)*((2*n+a+b)*(2*n+a+b-2*x + a^2 - b^2))
  part2 <- 2*(n+a-1)*(n+b-1)*(2*n+a+b)
  
  return((part1*JacobiPolynomial(n-1, a, b, x) - part2*JacobiPolynomial(n-2, a, b, x))/denom)
}

ChebyshevPolynomial <- function(n, x){
  if(n == 1) return(2*x)
  if(n == 0) return(1)
  
  return(2*x*ChebyshevPolynomial(n-1, x) - ChebyshevPolynomial(n-2, x))
}

