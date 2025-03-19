# Aux for uniform prediction band computation

library(doMC)
# source("CoDa_nonparametric_boot.R")
# source("CoDa_nonparametric_boot2.R")
packages <- c("generics", "demography", "forecast","fda","fdaoutlier", "rlist", "mrfDepth","ftsa","rainbow")
## Now load or install&load all
package_check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)



# Coda_nonparametric without differencing

### Eigenratio $k$ selection method ###
select_K <- function(tau, eigenvalue)
{
  
  k_max = length(eigenvalue)
  k_all = rep(0, k_max-1)
  for(k in 1:(k_max-1))
  {
    k_all[k] = (eigenvalue[k+1]/eigenvalue[k])*ifelse(eigenvalue[k]/eigenvalue[1] > tau, 1, 0) + ifelse(eigenvalue[k]/eigenvalue[1] < tau, 1, 0)
  }
  
  K_hat = which.min(k_all)
  
  return(K_hat)
}


# i=1
# dat=t(male_prefecture_res_means[[i]])
# fixed_comp=t(male_prefecture_fixed_means[[i]])
# # Point_forecasts_FM <- forecasted_FM_lrc_ARIMA_USA_4[[1]][[3]][,,1]
# alpha_transf=all_alpha_male[,i]
# sample_number=n_year
# fh=10
# B=50
# level=95
# fmethod ="arima"
# no_core=detectCores()-2
# ncomp_selection = "EVR"


# dat: p by n data matrix
# fh: forecast horizon
# B: number of bootstrap samples
# level: nominal coverage probability, such as level = 80
# fmethod: univariate time-series forecasting method

GW_LRC_nonstationary_boot_uniform<- function(dat,fixed_comp, alpha_transf,fh, B, level, fmethod = c("ets", "arima"),
                                             ncomp_selection = c("EVR", "fixed"))
{
  n = ncol(dat)
  p = nrow(dat)
  
  # compute mean
  
  my = fixed_comp[,1]
  mdata2 = array(rep(as.matrix(my), B * fh), dim = c(p, B, fh))
  
  # de-center data by subtracting the mean term
  
  new_dat = t(dat)
  
  LRC_est = long_run_covariance_estimation(dat = t(new_dat))
  eigen_decomp = eigen(LRC_est, symmetric = TRUE)
  
  # determine ncomp
  
  if(ncomp_selection == "EVR"){
    lambda_val = eigen_decomp$values[which(eigen_decomp$values > 0)]
    ncomp = select_K(tau = 1/log(max(lambda_val[1], length(lambda_val))), eigenvalue = lambda_val)
    rm(lambda_val); rm(LRC_est)
  }else{
    ncomp = 6
  }
  
  var_total_variations1 = (sum(eigen_decomp$values[1:ncomp])/sum(eigen_decomp$values))*100
  # compute the basis function and their scores
  
  LRC_basis = matrix(eigen_decomp$vectors[,1:ncomp], ncol = ncomp)
  LRC_score = new_dat %*% LRC_basis
  
  # obtain h-step-ahead forecast of insample principal component scores
  
  olivia = matrix(NA, ncomp, fh)
  if(fmethod == "ets")
  {
    for(i in 1:ncomp)
    {
      olivia[i,] = forecast(ets(LRC_score[,i]), h = fh)$mean
    }
  }else
  {
    for(i in 1:ncomp)
    {
      olivia[i,] = forecast(auto.arima(LRC_score[,i]), h = fh)$mean
    }
  }
  
  # compute insample forecast errors
  
  forerr = matrix(NA, (n - ncomp - fh + 1), ncomp)
  for(i in fh:(n - ncomp))
  {
    k = i + (ncomp - fh)
    fore = matrix(NA, 1, ncomp)
    if(fmethod == "ets")
    {
      for(j in 1:ncomp)
      {
        fore[,j] = forecast(ets(LRC_score[1:k, j]), h = fh)$mean[fh]
      }
    }
    else
    {
      for(j in 1:ncomp)
      {
        fore[,j] = forecast(auto.arima(LRC_score[1:k, j]), h = fh)$mean[fh]
      }
    }
    forerr[i - fh + 1,] = LRC_score[k + fh,] - fore
  }
  
  # compute functional residuals
  
  LRC_recon = LRC_basis %*% t(LRC_score)
  LRC_resi = t(new_dat) - LRC_recon
  
  # bootstrap the noise term
  
  q = array(NA, dim = c(p, B, fh))
  for(j in 1:fh)
  {
    for(i in 1:p)
    {
      q[i,,j] = sample(LRC_resi[i,], size = B, replace = TRUE)
    }
  }
  
  # bootstrap forecast errors of principal component scores
  
  ny = array(NA, dim = c(ncomp, B, fh))
  for(j in 1:fh)
  {
    for(i in 1:ncomp)
    {
      ny[i,,j] = sample(forerr[,i], size = B, replace = TRUE)
    }
  }
  
  # add the bootstrapped forecast errors to the forecast principal component scores
  
  oli = array(rep(olivia, B * fh), dim = c(ncomp, B, fh))
  fo = array(NA, dim = c(ncomp, B, fh))
  for(j in 1:fh)
  {
    for(i in 1:B)
    {
      fo[,i,j] = oli[,i,j] + ny[,i,j]
    }
  }
  
  # conditional on the estimated basis functions and mean, obtain forecast curves
  
  pred = array(NA, dim = c(p, B, fh))
  for(j in 1:fh)
  {
    for(i in 1:B)
    {
      pred[,i,j] = LRC_basis %*% fo[,i,j] + q[,i,j] + mdata2[,i,j]
    }
  }
  
  # commHS (please check the dimensionality)
  
  point_forecast = LRC_basis %*% olivia + mdata2[,1,]
  
  #transform back to the constrained space
  alpha_x_1  <- alpha_transf
  recon_1    <- (fixed_comp + dat)
  
  # # reconstruction (model in-sample fitting)
  f_x_t_star_recon_1 = d_x_t_star_recon_1 = matrix(NA, p, n)
  for(ik in 1:n)
  {
    f_x_t_star_recon_1[,ik] = exp(recon_1[,ik])/sum(exp(recon_1[,ik]))
    d_x_t_star_recon_1[,ik] = (f_x_t_star_recon_1[,ik] * alpha_x_1)/sum(f_x_t_star_recon_1[,ik] * alpha_x_1) * (10^5)
  }
  
  # back-transformation for the bootstrapped predictions
  boot_prev = boot_pred = array(NA, dim = c(p, B, fh))
  for(j in 1:fh)
  {
    for(i in 1:B)
    {
      boot_prev[,i,j] = exp(pred[,i,j])/sum(exp(pred[,i,j]))
      boot_pred[,i,j] = (boot_prev[,i,j] * alpha_x_1)/sum((boot_prev[,i,j] * alpha_x_1))*(10^5)
    }
  }
  
  # transform back the point forecasts
  PF_prev = matrix(NA, p,fh)
  PF_pred = matrix(NA, p,fh)
  for(j in 1:fh)
  {
      PF_prev[,j] = exp(point_forecast[,j])/sum(exp(point_forecast[,j]))
      PF_pred[,j] = (PF_prev[,j] * alpha_x_1)/sum((PF_prev[,j] * alpha_x_1))*(10^5)
  }
  
  # take corresponding quantiles
  
  k1 = k2 = k3 = matrix(NA, p, fh)
  for(j in 1:fh)
  {
    for(i in 1:p)
    {
      k1[i,j] = quantile(boot_pred[i,,j], (100 - level)/200, na.rm = TRUE)
      k2[i,j] = quantile(boot_pred[i,,j], 1 - (100 - level)/200, na.rm = TRUE)
      
      # commHS
      
      k3[i,j] = quantile((boot_pred[i,,j] - PF_pred[i,j]), level/100, na.rm = TRUE)
    }
  }
  
  # commHS
  
  uniform_band = cbind(PF_pred - k3, PF_pred + k3)
  
  lb <- PF_pred - k3
  ub <- PF_pred + k3
  
  # colnames(k1) <-c("h=1", "h=2", "h=3", "h=4", "h=5", "h=6", "h=7", "h=8", "h=9", "h=10")
  # colnames(k2) <-c("h=1", "h=2", "h=3", "h=4", "h=5", "h=6", "h=7", "h=8", "h=9", "h=10")
  
  return(list(bootsamp = pred, data_untransformed = d_x_t_star_recon_1, lb = lb, ub = ub,
              ncomp1 = ncomp, TV1 = var_total_variations1))
}

#############################################
# Function for computing the pointwise score
#############################################

# l: lower bound
# u: upper bound
# x: actual holdout data
# alpha: level of significance alpha = 0.2



build.tensor <- function(l, nrows, ncols,name,fh) {
  tensor <- array(0, dim = c(nrows, ncols, ncols))
  tensor[, , 1] <- l[[1]][[name]]
  for (i in 1:(fh-1)) {
    x <- fh - i + 1
    tensor[, -(x:fh), i+1] <- l[[i+1]][[name]]
  }
  tensor
}

compare.data_uniform <- function(lb=lb, ub=ub, real=test_data, fh=fh, level) {
  uniform_coverage = matrix(0, nrow = fh, ncol = 1)
  CPD_uniform <- matrix(0, nrow = fh, ncol = 1)
  
  for (i in 1:fh) {
    if (i == 1) {
      lb.point <- lb[, i, ]
      ub.point <- ub[, i, ]
    } else {
      x <- fh - i + 2
      lb.point <- as.matrix(lb[, i, -(x:fh)])
      ub.point <- as.matrix(ub[, i, -(x:fh)])
    }
    real.point <- as.matrix(real[, i:fh])
    lb_ind_count=ub_ind_count=count <-matrix(NA, nrow = dim(real.point)[1], ncol = dim(real.point)[2])
    for (j in 1:dim(real.point)[2]) {
      count[,j]=ifelse(any(c(any(which(lb.point[,j] >= real.point[,j])),
                   any(which(ub.point[,j] <= real.point[,j])))), 1, 0)
       # lb_ind_count[,j] = ifelse(any(real.point[,j] < lb.point[,j]), 1, 0)
       # ub_ind_count[,j] = ifelse(any(real.point[,j] > ub.point[,j]), 1, 0)
    }
    # lb_ind = ifelse(real.point < lb.point, 1, 0)
    # ub_ind = ifelse(real.point > ub.point, 1, 0)
    #Prod <-lb_ind*ub_ind
    # 
    uniform_coverage[i,] = 1 - sum(apply(count, 2, max))/ncol(real.point)
    CPD_uniform[i,]=abs(uniform_coverage[i,]-level/100) 
  }
  
  result = as.matrix(cbind(uniform_coverage,CPD_uniform))
  colnames(result) = c(paste("Uniform","coverage",level),paste("CPD Uniform", level, "% coverage"))
  return(result)
}



# 
# i=1
# dat=t(male_prefecture_res_means[[i]])
# fixed_comp=t(male_prefecture_fixed_means[[i]])
# alpha_transf=all_alpha_male[,i]
# sample_number=n_year
# fh=10
# B=50
# level=95
# fmethod ="arima"
# no_core=detectCores()-2
# ncomp_selection = "fixed"


coverage_CoDA_uniform <- function(dat,fixed_comp,alpha_transf,fh, B, level, fmethod ="arima",
                           sample_number=n_year,no_core,ncomp_selection = c("EVR", "fixed")){
  sim_data=dat
  grid_point = seq(0, 1, length.out = nrow(sim_data))
  colnames(sim_data) = 1:sample_number
  rownames(sim_data) = grid_point
  
  
  
  # define testing and training samples
  n_val = ncol(sim_data)
  n_training_ini = ncol(sim_data)-fh
  n_testing = n_val-n_training_ini
  
  # Construct the testing data
  d_x_t_star_recon_1=GW_LRC_nonstationary_boot_uniform(dat,fixed_comp, alpha_transf,fh=fh, B, level, fmethod,ncomp_selection)$data_untransformed
  test_data = as.matrix(d_x_t_star_recon_1[,(n_training_ini + 1):n_val])
  
  if((n_training_ini + n_testing) != n_val)
  {
    warning("length of training sample + testing sample != total sample")
  }else {
    
    # Compute the lb and ub at each forecast horizon.
    
    registerDoMC(no_core)
    boot = foreach(iwk = 1:n_testing) %dopar% GW_LRC_nonstationary_boot_uniform(dat=dat[,iwk:(n_training_ini - 1 + iwk)]
                                                                         ,fixed_comp=fixed_comp[,iwk:(n_training_ini - 1+ iwk)], 
                                                                         alpha_transf=alpha_transf,fh=(fh-iwk+1), B=B, level=level, 
                                                                         fmethod=fmethod,ncomp_selection =ncomp_selection)
    
    lb <- build.tensor(l=boot,nrows=nrow(sim_data),ncols=fh,name="lb",fh)
    ub <- build.tensor(l=boot,nrows=nrow(sim_data),ncols=fh,name="ub",fh)
    
    
    
  }
  
  result <-compare.data_uniform(lb=lb, ub=ub, real=test_data, fh=fh, level)
  
  return(list(lb=lb,ub=ub,result=result))
  
}





