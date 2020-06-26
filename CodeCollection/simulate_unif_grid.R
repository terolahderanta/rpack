
simulate_unif_grid = function(n_out, alpha=0.5, out_scale = 5, x=NULL, y=NULL){
  
  if(!is.null(x) & !is.null(y)){
  
    #abx = c(min(x), max(x))
    
    #aby = c(min(y), max(y))
    
    abx = c(min(x) - out_scale*sd(x)/sqrt(length(x)), max(x) + out_scale*sd(x)/sqrt(length(x)))
    
    aby = c(min(y) - out_scale*sd(y)/sqrt(length(y)), max(y) + out_scale*sd(y)/sqrt(length(y)))
    
  }
  
  i = 1
  
  d = 1
  
  while(d < n_out){
    
    x = runif(i*n_out, -1, 1)
    
    y = runif(i*n_out, -1, 1)
    
    x1 = x > alpha*max(x)
    
    x2 = x < alpha*min(x)
    
    y1 = y > alpha*max(y)
    
    y2 = y < alpha*min(y)
    
    X = x[x1 + x2 + y1 + y2 >= 1]
    
    Y = y[x1 + x2 + y1 + y2 >= 1]
    
    d = length(X)
    
    if(d > n_out){
      
      ind = sample(1:length(X), n_out)
      
      X = X[ind]
      
      Y = Y[ind]
      
    }
    
    i = i + 1
    
  }
  
  if(!is.null(x) & !is.null(y)){
    
    X = (abx[2] - abx[1])*((X - min(X))/(max(X) - min(X))) + abx[1]
    
    Y = (aby[2] - aby[1])*((Y - min(Y))/(max(Y) - min(Y))) + aby[1]
    
  }
  
  unif_grid_data = data.frame(x=X, y=Y)
  
  return(unif_grid_data)
  
}