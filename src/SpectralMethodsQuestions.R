source("include/ExponentialExample.R")
source("include/Quadrature.R")
source("include/DiscreteProjections.R")

t <- 5
f <- function(x) u(t, x)

x <- seq(-1, 1, length=100)
plot(x, f(x), type="l")
# for(i in 0:5){
#   lines(x, LegendreUniform(x, 20, f, N=i))
# }

lines(x, LegendreUniform(x, 20, f, N=0), col=2)
lines(x, LegendreUniform(x, 20, f, N=1), col=3)
lines(x, LegendreUniform(x, 20, f, N=2), col=4)
lines(x, LegendreUniform(x, 20, f, N=3), col=5)
lines(x, LegendreUniform(x, 20, f, N=4), col=6)
lines(x, LegendreUniform(x, 20, f, N=5), col=7)














































