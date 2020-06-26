random_sigma <- function (no_corr = TRUE){
  
  if (no_corr) {
    
    r_diag <- runif(1, min = 1, max = 1.5)
    
    Sigma <- diag(x = rep(r_diag, 2))
    
  }else{
    
    A <- matrix(runif(4, min = 0.1, max = 0.8), ncol = 2, 
                nrow = 2)
    
    Sigma <- t(A) %*% A
    
    Sigma[1, 1] <- Sigma[1, 1] + runif(1)
    
    Sigma[2, 2] <- Sigma[2, 2] + runif(1)
    
  }
  
  return(Sigma)
  
}