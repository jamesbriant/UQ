#######################################################################
# Exercise 1

library(mvtnorm)
#R package for simulating multivariate Gaussians

calcSigma <- function(X1,X2,theta){
  # Function takes in a discretised grid
  #input X1 x X2 and computes covariance matrix using covariance function k
  v1=theta[1]
  v2=theta[2]
  alpha=theta[3]
  r=theta[4]
  
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)){
    for (j in 1:ncol(Sigma)){
      Sigma[i,j] <- v1*exp(-(abs(X1[i]-X2[j])/r)^alpha) + v2*(i==j)
    }
  }
  
  return(Sigma)
}


theta=c(1,0.15,2,2)
f=data.frame(x=c(0.9,3.8,5.2,6.1,7.5,9.6), y=c(0.1,1.2,2.1,1.1,1.5,1.2))

par(mfrow=c(1,4))
plot(f$x,f$y,col=2,pch=10,ylim=c(-2.5,2.5),xlab="",ylab="")

x.star=seq(0,10,length=50)
Sig=calcSigma(x.star,x.star,theta)
prior.sample=matrix(0,nrow=length(x.star),ncol=20)

for(i in 1:ncol(prior.sample)){
  prior.sample[,i]=rmvnorm(1,mean=rep(0,length(x.star)), sigma=Sig)
}

matplot(x.star,prior.sample,type="l",xlab="",ylab="",main="prior sample",col="orange")

x=f$x
k.xx <- calcSigma(x,x,theta)
k.xxs <- calcSigma(x,x.star,theta)
k.xsx <- calcSigma(x.star,x,theta)
k.xsxs <- calcSigma(x.star,x.star,theta)

n.samples=20
# The standard deviation of the noise
sigma.n=theta[2]
# calculate the mean and covariance functions
f.star.bar <- k.xsx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%f$y
cov.f.star <- k.xsxs - k.xsx%*%solve(k.xx +sigma.n^2*diag(1, ncol(k.xx)))%*%k.xxs

# calculate the sample functions2
values <- matrix(0,nrow=length(x.star), ncol=n.samples)
for (i in 1:n.samples){
  values[,i] <- rmvnorm(1, f.star.bar, cov.f.star)
}

matplot(x.star,values,type="l",xlab="",ylab="",main="posterior sample",xaxs="i",yaxs="i",col="orange")
plot(x.star,f.star.bar,type="l",xlab="",ylab="",main="posterior mean+confidence band",col=4,ylim=c(-2.5,2.5))

y.up=f.star.bar+1.96*diag(cov.f.star)^0.5
y.low=f.star.bar-1.96*diag(cov.f.star)^0.5
points(f$x,f$y,col=2)
lines(x.star,y.low,col="dark grey",lty=2)
lines(x.star,y.up,col="dark grey",lty=2)





####### Prediction
x.predict <- seq(10, 11.5, length=10)

k.xpx <- calcSigma(x.predict,x,theta)
k.xxp <- calcSigma(x,x.predict,theta)
k.xpxp <- calcSigma(x.predict,x.predict,theta)

mu.star <- k.xpx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%f$y
sigma.star <- k.xpxp - k.xpx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%k.xxp

y.up=mu.star+1.96*diag(sigma.star)^0.5
y.low=mu.star-1.96*diag(sigma.star)^0.5

par(mfrow=c(1,1))
plot(x.predict, mu.star, ylim=c(min(y.low), max(y.up)))
points(x.predict, mu.star, type="l")
lines(x.predict, y.low,col="dark grey", lty=2)
lines(x.predict, y.up,col="dark grey", lty=2)




###########---------------------------------------------------------------
# Exercise 2
# Part 1

theta=c(1, 2, 0.03)
f=data.frame(x=c(0.9,3.8,5.2,6.1,7.5,9.6), y=c(0.1,1.2,2.1,1.1,1.5,1.2))

calcSigma <- function(X1,X2,theta){
  # Function takes in a discretised grid
  #input X1 x X2 and computes covariance matrix using covariance function k
  v1=theta[1]
  v2=theta[3]
  alpha=2
  r=theta[2]
  
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)){
    for (j in 1:ncol(Sigma)){
      Sigma[i,j] <- v1*exp(-(abs(X1[i]-X2[j])/r)^alpha) + v2*(i==j)
    }
  }
  
  return(Sigma)
}

x.plot <- seq(0, 10, length=50)

k.xx <- calcSigma(f$x, f$x, theta)
k.xplotx <- calcSigma(x.plot, f$x, theta)
k.xxplot <- calcSigma(f$x, x.plot, theta)
k.xplotxplot <- calcSigma(x.plot, x.plot, theta)

m.star <- k.xplotx %*% solve(k.xx + theta[3]*diag(1, ncol(k.xx))) %*% f$y
k.star <- k.xplotxplot - k.xplotx %*% solve(k.xx + theta[3]*diag(1, ncol(k.xx))) %*% k.xxplot

plot(x.plot, m.star, ylim=c(-1.25, 2.5))
lines(x.plot, m.star)
points(f$x, f$y, pch=16, col="red")

y.up <- m.star + 1.96*diag(k.star)^0.5
y.low <- m.star - 1.96*diag(k.star)^0.5
lines(x.plot, y.low, col="dark grey", lty=2)
lines(x.plot, y.up, col="dark grey", lty=2)


##########
# Part 2

MarginalLogLikelihood <- function(theta){
  f=data.frame(x=c(0.9,3.8,5.2,6.1,7.5,9.6), y=c(0.1,1.2,2.1,1.1,1.5,1.2))
  
  k.xx <- calcSigma(f$x, f$x, theta)
  
  a <- f$y %*% solve(k.xx) %*% f$y
  b <- determinant(k.xx, logarithm = TRUE)$modulus[1]
  c <- length(f$x) * log(2*pi)
  
  return(0.5*(a+b+c)) # optim performs minimisation by default, so negative the likelihood
}

theta.optimised <- optim(theta, MarginalLogLikelihood)
theta.optimised$par

v1.optimised=theta.optimised$par[1]
v2.optimised=theta.optimised$par[3]
r.optimised=theta.optimised$par[2]


#########
# Part 3

x.plot <- seq(0, 10, length=50)

k.xx.opt <- calcSigma(f$x, f$x, theta.optimised$par)
k.xplotx.opt <- calcSigma(x.plot, f$x, theta.optimised$par)
k.xxplot.opt <- calcSigma(f$x, x.plot, theta.optimised$par)
k.xplotxplot.opt <- calcSigma(x.plot, x.plot, theta.optimised$par)

m.star.opt <- k.xplotx.opt %*% solve(k.xx.opt + v2.optimised*diag(1, ncol(k.xx.opt))) %*% f$y
k.star.opt <- k.xplotxplot.opt - k.xplotx.opt %*% solve(k.xx.opt + v2.optimised*diag(1, ncol(k.xx.opt))) %*% k.xxplot.opt

points(x.plot, m.star.opt, col="blue")
lines(x.plot, m.star.opt, col="blue")
points(f$x, f$y, pch=16, col="red")

y.up.opt <- m.star.opt + 1.96*diag(k.star.opt)^0.5
y.low.opt <- m.star.opt - 1.96*diag(k.star.opt)^0.5
lines(x.plot, y.low.opt, col="blue", lty=2)
lines(x.plot, y.up.opt, col="blue", lty=2)


#########---------------------------------------------------------
# Exercise 3

library(deSolve)

alpha <- 1
beta <- 2
gamma <- 3
delta <- 4

x0 <- 2
y0 <- 3

tt <- (1:100)/10

intial.conditions <- c(x0, y0)

odes <- function(t, data, params){
  alpha <- params[1]
  beta <- params[2]
  gamma <- params[3]
  delta <- params[4]
  
  x <- data[1]
  y <- data[2]
  
  xprime <- alpha*x - beta*x*y
  yprime <- delta*x*y - gamma*y
  
  return(list(c(xprime, yprime)))
}

solution <- as.data.frame(ode(intial.conditions, tt, odes, c(alpha, beta, gamma, delta)))
names(solution) <- c("time", "x", "y")
plot(tt, solution$x)
points(tt, solution$y, col="red")

# Add noise to the data
D <- solution
D$x <- D$x + rnorm(length(tt), 0, 0.3)
D$y <- D$y + rnorm(length(tt), 0, 0.5)



