library(deSolve)

logistic=function(t, y, params){
  c = params[1]
  r = params[2]
  dy = r*y*(1-y/c)
  
  return(list(dy))
}

tt <- 1:10
y_init <- 2
init.params = c(c=100, r=2)#  c=100,r=2

## use function `rk' to solve ODE with Runge-Kutta method
## output contains days in 1st col, growth in 2nd col

sim.data = as.data.frame(rk(y_init, tt, logistic, init.params))

## Add noise to solution to generate data
set.seed(12345)

sim.data$noisy = sim.data[, 2] + rnorm(nrow(sim.data), 0, 13)
colnames(sim.data) = c("day", "growth", "noisy.growth")

## Now carry out nonlinear regression with data with the logistic model
param.start = list(c = max(sim.data$noisy.growth), r=2)
nls.fit = nls(noisy.growth ~ rk(noisy.growth, day, logistic, c(c=c, r=r))[, 2], data=sim.data, start=param.start)

### To get estimate of c and r parameters, their variances, and fitted curve
coef(nls.fit)

##         c         r
## 96.311270  1.265456

vcov(nls.fit)

##            c           r
## c 20.6717037 -0.31755687
## r -0.3175569  0.02772272

nls.estimate = fitted(nls.fit)

plot(noisy.growth ~ day, data=sim.data, xlab="day", ylab="growth")
lines(noisy.growth ~ day, data=sim.data)
lines(growth ~ day, data=sim.data, col="blue")
lines(sim.data$day, nls.estimate, col="red")
legend('bottomright', legend = c("data", "ODE solution", "NLS fitted curve"), col = c("black", "blue", "red"), pt.lwd=1, lwd=1, bty="n", text.font=.5,border=F)



