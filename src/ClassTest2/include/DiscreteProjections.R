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





BetaJacobi <- function(lambda, J, user.func, N=0, lambda.dagger=0, epsilon=1){
  # lambda.dagger and epsilon map [-1, 1] on to the region of interest
  #   lambda: evaluation points, (these are pre-determined, NOT random)
  #   J: Number of quadrature points for Composite Simpson's Rule
  #   user.func: Must accept lambda only (and not time)
  #   N: Polynomial size
  
  
  approximation <- numeric(length(lambda))
  u.n <- numeric(N+1)
  gammas <- numeric(N+1)
  for(n in 0:N){
    # evaluate the gamma
    gamma <- GetGammas(user.func=function(x, i=n) JacobiPolynomial(i, 1.5, 1.5, x),
                       dens.func=function(x) return(dbeta(x, 1.5, 1.5)),
                       a=0, b=1)
    
    gammas[n+1] <- gamma
    
    # calculate the u term
    u <- CompositeSimpsonRule(function(x) user.func(x)*JacobiPolynomial(i, 1.5, 1.5, x), 
                              N=J, a=-1, b=1)/gamma
    
    # save the coefficients
    u.n[n+1] <- u
    
    # 
    approximation <- approximation + u*JacobiPolynomial(n, 1.5, 1.5, lambda)
  }
  
  
  return(list(approximation=approximation,
              coefficients=u.n,
              gammas=gammas)
  )
  
}













WignerChebyshev <- function(lambda, J, user.func, N=0){
  #   lambda: evaluation points, (these are pre-determined, NOT random)
  #   J: Number of quadrature points for Composite Simpson's Rule
  #   user.func: Must accept lambda only (and not time)
  #   N: Polynomial size
  
  approximation <- numeric(length(lambda))
  u.n <- numeric(N+1)
  gammas <- numeric(N+1)
  for(n in 0:N){
    # evaluate the gamma
    gamma <- GetGammas(user.func=function(x, i=n) ChebyshevPolynomial(i, x),
                       dens.func=dwigner,
                       a=-1, b=1)*2
    gammas[n+1] <- gamma
    
    # calculate the u term
    u <- CompositeSimpsonRule(function(x) user.func(x)*ChebyshevPolynomial(n, x), 
                              N=J, a=-1, b=1)/gamma
    
    # save the coefficients
    u.n[n+1] <- u
    
    # 
    approximation <- approximation + u*ChebyshevPolynomial(n, lambda)
  }
  
  
  return(list(approximation=approximation,
              coefficients=u.n,
              gammas=gammas)
         )
}