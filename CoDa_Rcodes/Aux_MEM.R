#Auxiliary functions for computing the forecast errors in the MEM method. 

fh=10
extract_forecast <- function(i, h, tensor) {
  tensor[, i, h-i+1] 
}

source("forecast_errors.R")
compute_error <- function(forecast, true){
  
  
  MAPE   <- mape(forecast, true)
  KLD    <- KLD(forecast, true)
  JS1    <- JSD_simple(forecast,true)
  JS2    <- JSD_geom(forecast, true)
  
  
  return(c(MAPE,KLD,JS1,JS2))
}


pairwise_errors <- function(j, forecast, true) {
  forecast <- forecast[, j]
  if (is.null(dim(true))) true <- matrix(true, ncol = 1)
  true <- true[, j]
  compute_error(forecast, true)
}

build.tensor <- function(l, nrows, ncols,name) {
  tensor <- array(0, dim = c(nrows, ncols, ncols))
  tensor[, , 1] <- l[[1]][[name]]
  for (i in 1:(ncols-1)) {
    x <- ncols - i + 1
    tensor[, -(x:ncols), i+1] <- l[[i+1]][[name]]
  }
  tensor
}
