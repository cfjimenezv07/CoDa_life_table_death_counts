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
# dat=t(male_prefecture_res[[i]])
# fixed_comp=t(male_prefecture_fixed[[i]])
# alpha_transf=all_alpha_male[,i]
# sample_number=n_year
# fh=10
# B=10
# level=80
# fmethod ="arima"
# no_core=detectCores()-2
# ncomp_selection = "EVR"


# dat: p by n data matrix
# fh: forecast horizon
# B: number of bootstrap samples
# level: nominal coverage probability, such as level = 80
# fmethod: univariate time-series forecasting method

GW_LRC_nonstationary_boot2 <- function(dat,fixed_comp, alpha_transf,fh, B, level, fmethod = c("ets", "arima"),ncomp_selection = c("EVR", "fixed"))
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
  
  # obtain h-step-ahead forecast of principal component scores
  
  olivia = matrix(NA, ncomp, fh)
  if(fmethod == "ets")
  {
    for(i in 1:ncomp)
    {
      olivia[i,] = forecast(ets(LRC_score[,i]), h = fh)$mean
    }
  }
  else
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
    else{
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
  
  
  
  #transform back to the constrained space
  alpha_x_1  <- alpha_transf
  recon_1    <- (fixed_comp+dat)
  
  # # reconstruction (model in-sample fitting)
  f_x_t_star_recon_1 = d_x_t_star_recon_1 = matrix(NA, p, n)
  for(ik in 1:n)
  {
    f_x_t_star_recon_1[,ik] = exp(recon_1[,ik])/sum(exp(recon_1[,ik]))
    d_x_t_star_recon_1[,ik] = (f_x_t_star_recon_1[,ik] * alpha_x_1)/sum(f_x_t_star_recon_1[,ik] * alpha_x_1)*(10^5)
  }
  
  
  # back-transformation for the bootstrapped predictions
  boot_prev <- array(NA, dim = c(p, B, fh))
  boot_pred <- array(NA, dim = c(p, B, fh))
  for(j in 1:fh)
  {
    for(i in 1:B)
    {
      boot_prev[,i,j] = exp(pred[,i,j])/sum(exp(pred[,i,j]))
      boot_pred[,i,j] = (boot_prev[,i,j] * alpha_x_1)/sum((boot_prev[,i,j] * alpha_x_1))*(10^5)
    }
  }
  
  
  # take corresponding quantiles
  
  k1 = k2 = matrix(NA, p, fh)
  for(j in 1:fh)
  {
    for(i in 1:p)
    {
      k1[i,j] = quantile(boot_pred[i,,j], (100 - level)/200, na.rm = TRUE)
      k2[i,j] = quantile(boot_pred[i,,j], 1 - (100 - level)/200, na.rm = TRUE)
    }
  }
  # colnames(k1) <-c("h=1","h=2",'h=3',"h=4","h=5","h=6","h=7","h=8","h=9","h=10")
  # colnames(k2) <-c("h=1","h=2",'h=3',"h=4","h=5","h=6","h=7","h=8","h=9","h=10")
  
  return(list(bootsamp = pred,data_untransformed=d_x_t_star_recon_1, lb = k1, ub = k2,ncomp1=ncomp,TV1=var_total_variations1))
}
