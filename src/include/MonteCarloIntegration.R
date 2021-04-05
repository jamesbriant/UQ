MonteCarlo2D <- function(user.func, sample=FALSE, a0=0, a1=1, b0=0, b1=1){
  # Assumes domain is [0,1]x[0,1]
  # domain is rectangular
  
  if(identical(sample, FALSE)){
    sample = Generate2DUniformSample(N=5000, a0, a1, b0, b1)
  }
  
  data <- user.func(sample$x, sample$y)
  
  return(list(mean=mean(data),
              variance=var(data)
              )
         )
}

Generate2DUniformSample <- function(N=1000, a0=0, a1=1, b0=0, b1=1){
  return(list(x=runif(N, a0, a1), 
              y=runif(N, b0, b1)
              )
         )
}





