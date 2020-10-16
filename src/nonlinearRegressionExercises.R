library(deSolve)

f <- function(x, c, r){
  a <- (c-2)/2
  
  return(c/(1+a*exp(-r*x)))
}

logistic=function(t, y, params){
  c = params[1]
  r = params[2]
  dy = r*y*(1-y/c)
  
  return(list(dy))
}

c.true <- 100
r.true <- 2

x <- 1:10
fx <- f(x, c.true, r.true)
set.seed(161020)
fx.noisy <- fx + rnorm(length(x), 0, 5^2)

nls.fit <- nls(fx.noisy ~ f(x, c, r), start=list(c=95, r=1.8))



coef(nls.fit)
vcov(nls.fit)

nls.estimate = fitted(nls.fit)

plot(x, fx.noisy, xlab="day", ylab="growth")
lines(x, fx.noisy)
lines(x, fx, col="blue")
lines(x, nls.estimate, col="red")
legend('bottomright', legend = c("data", "true solution to ODE", "NLS fitted curve"), col = c("black", "blue", "red"), pt.lwd=1, lwd=1, bty="n", text.font=.5,border=F)





