# auxiliary functions for the TNH method

##### set a working directory
setwd("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@kaust.edu.sa/My Drive/Spring_2024/STAT_397/PhD_Project_5/Codes_paper2")
source('./MITS_class.R') 

init_data_new <-  function(data,p,Names){
  NAMES = Names
  N <- length(data)
  out = c('data','mixed_fts')
  mixed_fts <- mits$new()
  
  age_max = dim(data[[1]])[1]
  T = dim(data[[1]])[2]
  basis <- create.bspline.basis(c(0, 1), nbasis = p, norder = 4)
  args <- seq(0, 1, length=age_max)
  
  fd<-list()
  for(i in 1:N){
    fd[[i]] <- Data2fd(args, data[[i]], basis)
    mixed_fts$loadFunc(fd[[i]])
  }
  mixed_fts$actualize()
  names(fd) <- NAMES
  names(data) <- NAMES
  return(mget(out))
}

trim_data_new <- function(data_frame,T){
  out <- c('data')
  data <- list()
  names <- names(data_frame)
  for(i in 1:length(data_frame)){
    data[[names[i]]] <- data_frame[[names[i]]][,T:(T+n_year-h_max-1)]
    rownames(data[[names[i]]]) <- rownames(data_frame[[names[i]]])
    colnames(data[[names[i]]]) <- colnames(data_frame[[names[i]]])[T:(T+n_year-h_max-1)]
  }
  return(mget(out))
}


# Remove the zeros
remove_zeroes <- function(data_raw,n_prefectures,n_year,n_age) {
  N       <- n_prefectures
  T       <- n_year
  age_max <- n_age
  
  data_nozero <- list()
  for(i in 1:N) {
    data_nozero[[i]] <- matrix(data_raw[[i]],age_max,T)
  }
  for(i in 1:N) {
    for(j in 2:age_max) {#For j=1, the two zeroes have been removed before
      for(k in 1:T){ 
        if(data_nozero[[i]][j,k]==-Inf){
          data_nozero[[i]][j,k] <- data_nozero[[i]][j-1,k]
        }
      }
    }
  }
  return(data_nozero)
}

transform_back <- function(j,tensor, max_h,alpha_x,n_year,n_age,age){
  
  fore_val        <-  tensor[,,j]
  # back-transformation
  
  f_x_t_star_fore = d_x_t_star_fore = matrix(NA, n_age, max_h)
  for(ik in 1:max_h)
  {
    f_x_t_star_fore[,ik] = exp(fore_val[,ik])/sum(exp(fore_val[,ik]))
    d_x_t_star_fore[,ik] = (f_x_t_star_fore[,ik] * alpha_x)/sum((f_x_t_star_fore[,ik] * alpha_x))
  }
  colnames(d_x_t_star_fore)  = 1:max_h
  rownames(d_x_t_star_fore)  = age
  # 
  Forc_transformed <- d_x_t_star_fore*10^5; if (j > 1) Forc_transformed[,-((max_h-j+2):max_h)] -> Forc_transformed
  
  return(list(Forc_transformed = Forc_transformed))
  
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


extract_forecast <- function(i, h, tensor) {
  tensor[, i, h-i+1] 
}

source("forecast_errors.R")
compute_error <- function(forecast, true){
  
  
  KLD    <- KLD(forecast, true)
  JS1    <- JSD_simple(forecast,true)
  JS2    <- JSD_geom(forecast, true)
  
  
  return(c(KLD,JS1,JS2))
}


pairwise_errors <- function(j, forecast, true) {
  forecast <- forecast[, j]
  if (is.null(dim(true))) true <- matrix(true, ncol = 1)
  true <- true[, j]
  compute_error(forecast, true)
}

build.tensor2 <- function(l, nrows, ncols) {
  tensor <- array(0, dim = c(nrows, ncols, ncols))
  tensor[, , 1] <- l[[1]]
  for (i in 1:(ncols-1)) {
    x <- ncols - i + 1
    tensor[, -(x:ncols), i+1] <- l[[i+1]]
  }
  tensor
}

norm_col <- function(data){
  n_age <- dim(data)[1]
  h_max <- dim(data)[2]
  new_data <- matrix(0,n_age,h_max)
    for (j in 1:h_max) {
        new_data[,j] <-(data[,j]/sum(data[,j]))*(10^5)
      }

  return(new_data)
}


Remove_zeros_norm <-function(data){
  n_age <- dim(data)[1]
  h_max <- dim(data)[2]
  new_data <- matrix(0,n_age,h_max)
  for (i in 1:n_age) {
    for (j in 1:h_max) {
      if(data[i,j]<0){
        new_data[i,j] <- 0
      }else{
        new_data[i,j] <-data[i,j]
      }
    }
  }
  Norm_data <- norm_col(new_data)
  return(Norm_data)
}


