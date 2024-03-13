# Point forecast computation

#General function for point forecasts based on the FMP-ANOVA/FM-ANOVA decomposition

#1. General function for point forecasts obtained  to the functional residuals
# after removing deterministic components from the FMP-ANOVA/FM-ANOVA approach

#2. Computation of the point forecasts based on the rolling window approach, 
# using the function Pref_forecast_curves.


###############################################################################
# Computation of the point forecast based on the FMP decomposition
################################################################################
library("demography")
source("forecast_Arima.R")
library("ftsa")
library("vars")
library("doMC")
select_k <- function(tau, eigenvalue)
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

# #Example
# i=1
# pref=(n_year * i - (n_year - 1)):(n_year * i)
# k=10
# pref_k <- pref[k:(n_year - 10 - 1 + k)]
# 
# fixed_com = Fixed_part[pref_k, ]
# Residuals_f = Residuals_[pref_k, ]


Pref_forecast_curves<-function(fixed_com,Residuals_f,
                               est_method = c("lrc", "cov"),
                               fh = 30, 
                               B = 1000, 
                               prediction_method=c("ARIMA","VAR"),select_K=c("Fixed","EVR"), K=6){
  med_polish_resi=t(Residuals_f)
  if(est_method == "lrc"){
    # estimate long-run covariance by kernel sandwich estimator
    med_polish_resi_lrc = long_run_covariance_estimation(med_polish_resi)
  }else if(est_method == "cov"){
    # estimate empirical covariance function
    med_polish_resi_lrc = cov(t(med_polish_resi))
  }
  # perform eigen-decomposition
  med_polish_resi_eigen = eigen(med_polish_resi_lrc)
  
  if(select_K=="Fixed"){
    
    retain_component = K
    
  }else if(select_K=="EVR"){
    # determine retained number of components via eigenvalue ratio
    
    lambda_val = med_polish_resi_eigen$values[which(med_polish_resi_eigen$values > 0)]
    retain_component = select_k(tau = 1/log(max(lambda_val[1], length(lambda_val))), eigenvalue = lambda_val)
    
    
  }
  
  var_total_variations = (sum(med_polish_resi_eigen$values[1:retain_component])/sum(med_polish_resi_eigen$values))*100
  
  
  # determine 1st set of basis function and its scores
  med_polish_resi_basis = as.matrix(med_polish_resi_eigen$vectors[,1:retain_component])
  med_polish_resi_score = crossprod(med_polish_resi, med_polish_resi_basis)
  
  # obtain forecasts of PC scores via auto.arima
  med_polish_resi_score_forecast = matrix(NA, retain_component, fh)
  med_polish_resi_score_forecast_boot = array(NA, dim = c(retain_component, fh, B))
  if(prediction_method=="ARIMA"){
    for(ik in 1:retain_component){
      dum = forecast_Arima(object=auto.arima(med_polish_resi_score[,ik]), h = fh, bootstrap = TRUE, npaths = B)
      med_polish_resi_score_forecast[ik,] = dum$mean
      med_polish_resi_score_forecast_boot[ik,,] = t(dum$sim)
      # rm(ik); rm(dum)
    }
  }else if(prediction_method=="VAR"){
    object=med_polish_resi_score
    colnames(object)<-1:dim(object)[2]
    lag=VARselect(y=object,type = "const")$selection[1]
    model_VAR <- VAR(y=object,type = "const",ic="AIC",p=lag)
    pred=predict(model_VAR,n.ahead=fh)$fcst
    for (ik in 1:retain_component) {
      pred1=pred[[ik]]
      med_polish_resi_score_forecast[ik,]=pred1[,1]
    }
    
  }
  
  
  med_polish_resi_forecast = med_polish_resi_basis %*% med_polish_resi_score_forecast
  
  # add the fixed parts
  
  Fixed=t(fixed_com)[,1:fh]
  med_polish_curve_forecast = med_polish_resi_forecast + Fixed
  
  return(list(med_polish_curve_forecast=med_polish_curve_forecast, 
              med_polish_resi_forecast=  med_polish_resi_forecast,TV = var_total_variations))
  
}

build.tensor <- function(l, nrows, ncols,name) {
  tensor <- array(0, dim = c(nrows, ncols, ncols))
  tensor[, , 1] <- l[[1]][[name]]
  for (i in 1:(max_h-1)) {
    x <- max_h - i + 1
    tensor[, -(x:max_h), i+1] <- l[[i+1]][[name]]
  }
  tensor
}

################################################################################
# With rolling window approach
################################################################################
transform_back <- function(j,tensor, max_h,alpha_x_1,alpha_x_2,n_year,n_age,age){

  forecasted_curves        <-  tensor[,,j]
  forecasted_curves_male   <- forecasted_curves[1:n_age,]
  forecasted_curves_female <- forecasted_curves[(n_age+1):(2*n_age),]
  fore_val_1 <- forecasted_curves_male
  fore_val_2 <- forecasted_curves_female
  
  # back-transformation
  
  f_x_t_star_fore_1 = f_x_t_star_fore_2 = 
    d_x_t_star_fore_1 = d_x_t_star_fore_2 = matrix(NA, n_age, max_h)
  for(ik in 1:max_h)
  {
    f_x_t_star_fore_1[,ik] = exp(fore_val_1[,ik])/sum(exp(fore_val_1[,ik]))
    f_x_t_star_fore_2[,ik] = exp(fore_val_2[,ik])/sum(exp(fore_val_2[,ik]))
    d_x_t_star_fore_1[,ik] = (f_x_t_star_fore_1[,ik] * alpha_x_1)/sum((f_x_t_star_fore_1[,ik] * alpha_x_1))
    d_x_t_star_fore_2[,ik] = (f_x_t_star_fore_2[,ik] * alpha_x_2)/sum((f_x_t_star_fore_2[,ik] * alpha_x_2))
  }
  colnames(d_x_t_star_fore_1) = colnames(d_x_t_star_fore_2) = 1:max_h
  rownames(d_x_t_star_fore_1) = rownames(d_x_t_star_fore_2) = age
  # 
  Forc_transformed_male <- d_x_t_star_fore_1*10^5;   if (j > 1) Forc_transformed_male[,-((max_h-j+2):max_h)] -> Forc_transformed_male
  Forc_transformed_female <- d_x_t_star_fore_2*10^5; if (j > 1) Forc_transformed_female[,-((max_h-j+2):max_h)] -> Forc_transformed_female

  return(list(forecast_male = Forc_transformed_male,forecast_female = Forc_transformed_female))
  
}



# i=10
# fixed_com=Fixed_part
# Residuals_f=Residuals_
# est_method = "lrc"
# prediction_method="ARIMA"
# select_K="Fixed"
# K=6
# all_alpha_male=all_alpha_male
# all_alpha_female=all_alpha_female
# n_year=62
# n_age=111
# age=age
# no_core=detectCores()-2

library(foreach)
library(doParallel)
max_h = 10
ForecastC <- function(i,fixed_com,Residuals_f,est_method = c("lrc", "cov"),prediction_method=c("ARIMA","VAR"),
                      select_K=c("Fixed","EVR"), K=6,no_core=detectCores()-2,
                      all_alpha_male,all_alpha_female,n_year,n_age,age)
{
  pref=(n_year * i - (n_year - 1)):(n_year * i)
  n_training_ini =length(pref)-max_h
  fixed_com1 = fixed_com[pref,]
  Residuals_f1 = Residuals_f[pref,]
  registerDoMC(no_core)
  boot = foreach(iwk = 1:max_h ) %dopar% Pref_forecast_curves(fixed_com = fixed_com1[iwk:(n_training_ini - 1 + iwk),],
                                                                 Residuals_f = Residuals_f1[iwk:(n_training_ini - 1 + iwk),],
                                                                 est_method = est_method,
                                                                 fh = (max_h -iwk+1), B = 1000, 
                                                                 prediction_method=prediction_method,select_K=select_K, K=K)
  
  # Create the tensor
  Tensor <- build.tensor(l=boot, nrows=2*n_age, ncols=max_h,name="med_polish_curve_forecast")

  
  # all the data for each prefecture
  Y = fixed_com1 + Residuals_f1
  alpha_x_1  <- all_alpha_male[,i]
  alpha_x_2  <- all_alpha_female[,i]
  recon_1    <- t(Y[,1:n_age])
  recon_2    <- t(Y[,(n_age+1):(2*n_age)])
  # # reconstruction (model in-sample fitting)
  f_x_t_star_recon_1 = f_x_t_star_recon_2 =
    d_x_t_star_recon_1 = d_x_t_star_recon_2 = matrix(NA, n_age, n_year)
  for(ik in 1:n_year)
  {
    f_x_t_star_recon_1[,ik] = exp(recon_1[,ik])/sum(exp(recon_1[,ik]))
    f_x_t_star_recon_2[,ik] = exp(recon_2[,ik])/sum(exp(recon_2[,ik]))
    d_x_t_star_recon_1[,ik] = (f_x_t_star_recon_1[,ik] * alpha_x_1)/sum(f_x_t_star_recon_1[,ik] * alpha_x_1)
    d_x_t_star_recon_2[,ik] = (f_x_t_star_recon_2[,ik] * alpha_x_2)/sum(f_x_t_star_recon_2[,ik] * alpha_x_2)
  }
  
  dat_transformed_male   <- d_x_t_star_recon_1*10^5 
  dat_transformed_female <- d_x_t_star_recon_2*10^5
  
  # transform the tensor to the orignal space.
  
  Transform_tensor <-lapply(1:max_h, transform_back , Tensor, max_h,alpha_x_1,alpha_x_2,n_year,n_age,age)

  Forc_transformed_male   <- build.tensor(l=Transform_tensor, nrows=n_age, ncols=max_h,name="forecast_male")
  Forc_transformed_female <- build.tensor(l=Transform_tensor, nrows=n_age, ncols=max_h,name="forecast_female")
  
  print(i)
  return(list(dat_transformed_male =dat_transformed_male ,dat_transformed_female=dat_transformed_female,
              Forc_transformed_male=Forc_transformed_male,Forc_transformed_female=Forc_transformed_female))
}


