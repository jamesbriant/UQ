StepFunction2D <- function(x, y){
  # Domain is [0, 1]
  
  centre.x <- c(0.4, 0.6)
  centre.y <- c(0.4, 0.6)
  
  n <- length(x)
  output <- numeric(n)
  
  for(i in 1:n){
    if(x[i] > centre.x[1] && x[i] < centre.x[2]){
      if(y[i] > centre.y[1] && y[i] < centre.y[2]){
        output[i] <- 1
      }else{
        output[i] <- 0
      }
    }else{
      output[i] <- 0
    }
  }
  
  return(output)
}

