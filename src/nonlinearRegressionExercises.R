########------------------------------------------------------------------
# Exercise 1

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

nls.fit <- nls(fx.noisy ~ f(x, c, r), start=list(c=100, r=2))



coef(nls.fit)
vcov(nls.fit)

nls.estimate = fitted(nls.fit)

plot(x, fx.noisy, xlab="day", ylab="growth")
lines(x, fx.noisy)
lines(x, fx, col="blue")
lines(x, nls.estimate, col="red")
legend('bottomright', legend = c("data", "true solution to ODE", "NLS fitted curve"), col = c("black", "blue", "red"), pt.lwd=1, lwd=1, bty="n", text.font=.5,border=F)



#########--------------------------------------------------------------
# Exercise 2

library(deSolve)

f <- function(conc, theta1, theta2){
  return(theta1*conc/(theta2+conc))
}

conc <- c(2.856829, 5.005303, 7.519473, 22.101664, 27.769976, 39.198025, 45.483269, 203.784238)
rate <- c(14.58342, 24.74123, 31.34551, 72.96985, 77.50099, 96.08794, 96.96624, 108.88374)
data <- data.frame(conc, rate)

nls.fit <- nls(rate ~ f(conc, theta1, theta2), data=data, start=list(theta1=100, theta2=15))

nls.estimate = predict(nls.fit)

plot(conc, rate)
lines(conc, rate)
lines(conc, nls.estimate, col="red")
legend('bottomright', legend = c("data", "NLS fitted curve"), col = c("black", "red"), pt.lwd=1, lwd=1, bty="n", text.font=.5,border=F)











