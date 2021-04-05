GetGammas <- function(user.func, dens.func, a=0, b=1){
  CompositeSimpsonRule(function(x) return(user.func(x)^2*dens.func(x)),
                       N=20,
                       a=a, b=b
                       )
}

LegendreUniform <- function(lambda, J, user.func, N=0, lambda.dagger=0, epsilon=1){
  # lambda.dagger and epsilon map [-1, 1] on to the region of interest
  #   lambda: evaluation points, (these are pre-determined, NOT random)
  #   J: Number of quadrature points for Composite Simpson's Rule
  #   user.func: Must accept lambda only (and not time)
  #   N: Polynomial size
    
  
  approximation <- numeric(length(lambda))
  for(n in 0:N){
    # evaluate the gamma
    gamma <- GetGammas(user.func=function(x, i=n) LegendrePolynomial(i, x),
                       dens.func=function(x) return(0.5),
                       a=-1, b=1)
    
    # calculate the u term
    u <- CompositeSimpsonRule(function(x) user.func(x)*LegendrePolynomial(n, x), 
                              N=J, a=-1, b=1)/gamma
    
    # 
    approximation <- approximation + u*LegendrePolynomial(n, lambda)
  }
  

  return(approximation)
  
}