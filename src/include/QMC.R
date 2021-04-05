# Van der Corput implementation
VDC <- function(n, base=2){
  vdc <- 0
  denom <- 1
  
  while(n){
    denom <- denom*base
    remainder <- n %% base
    n <- n %/% base
    vdc <- vdc + remainder/denom
  }
  
  return(vdc)
}

# Returns the first N terms of the Van der Corput sequence base=base
VDCSequence <- function(N, base=2){
  sequence <- numeric(N)
  for(i in 1:N){
    sequence[i] <- VDC(i, base=base)
  }
  
  return(sequence)
}



GCD <- function(p, q){
  # Create the gcd of two positive integers.
  while(q != 0){
    temp <- p
    p <- q
    q <- temp %% q
  }
  
  return(p)
}

IsCoprime <- function(x, y) return(GCD(x, y) == 1)


# Returns the first N terms of the Halton sequence with bases base1 and base2
HaltonSequence <- function(N, base1=2, base2=3){
  if(!IsCoprime(base1, base2)){
    print("WARNING: base1 and base2 must be coprime.")
  }
  
  sequence.base1 <- VDCSequence(N, base=base1)
  sequence.base2 <- VDCSequence(N, base=base2)
  
  return(list("base1"=sequence.base1, 
              "base2"=sequence.base2)
         )
}




