library(doMC)
# source("CoDa_nonparametric_boot.R")
source("CoDa_nonparametric_boot2.R")
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
#############################################
# Function for computing the pointwise score
#############################################

# l: lower bound
# u: upper bound
# x: actual holdout data
# alpha: level of significance alpha = 0.2

interval_score <- function(holdout, lb, ub, alpha,na.rm=TRUE)
{
  lb_ind = ifelse(holdout < lb, 1, 0)
  ub_ind = ifelse(holdout > ub, 1, 0)
  score  = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
  return(score)
}

build.tensor <- function(l, nrows, ncols,name,fh) {
  tensor <- array(0, dim = c(nrows, ncols, ncols))
  tensor[, , 1] <- l[[1]][[name]]
  for (i in 1:(fh-1)) {
    x <- fh - i + 1
    tensor[, -(x:fh), i+1] <- l[[i+1]][[name]]
  }
  tensor
}

compare.data <- function(lb=lb, ub=ub, real=test_data, fh=fh, level) {
  pointwise_coverage <- matrix(0, nrow = fh, ncol = 1)
  pointwise_score <- matrix(0, nrow = fh, ncol = 1)
  CPD_pointwise <- matrix(0, nrow = fh, ncol = 1)
  
  for (i in 1:fh) {
    if (i == 1) {
      lb.point <- lb[, i, ]
      ub.point <- ub[, i, ]
    } else {
      x <- fh - i + 2
      lb.point <- lb[, i, -(x:fh)]
      ub.point <- ub[, i, -(x:fh)]
    }
    real.point <- real[, i:fh]
    lb_ind = which(lb.point >= real.point)
    ub_ind = which(ub.point <= real.point)
    #Prod <-lb_ind*ub_ind
    pointwise_coverage[i,] = 1 - length(union(lb_ind, ub_ind))/length(real.point)
    pointwise_score[i,] = mean(interval_score(holdout = real.point, lb = lb.point,
                                              ub = ub.point, alpha = (1-level/100)))
    
    CPD_pointwise[i,]=abs(pointwise_coverage[i,]-level/100) 
  }
  result = as.matrix(cbind(pointwise_coverage,CPD_pointwise,pointwise_score))
  colnames(result) = c(paste("Pointwise","coverage",level),paste("CPD Pointwise", level, "% coverage"),paste(level,"%", "interval score"))
  return(result)
}

# i=1
# dat=t(male_prefecture_res[[i]])
# fixed_comp=t(male_prefecture_fixed[[i]])
# alpha_transf=all_alpha_male[,i]
# sample_number=n_year
# fh=10
# B=10
# level=95
# fmethod ="arima"
# no_core=detectCores()-2
# ncomp_selection = "fixed"


coverage_CoDA2 <- function(dat,fixed_comp,alpha_transf,fh, B, level, fmethod ="arima",
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
  d_x_t_star_recon_1=GW_LRC_nonstationary_boot2(dat,fixed_comp, alpha_transf,fh=fh, B, level, fmethod,ncomp_selection)$data_untransformed
  test_data = as.matrix(d_x_t_star_recon_1[,(n_training_ini + 1):n_val])
  
  if((n_training_ini + n_testing) != n_val)
  {
    warning("length of training sample + testing sample != total sample")
  }else {
    
    # Compute the lb and ub at each forecast horizon.

    registerDoMC(no_core)
    boot = foreach(iwk = 1:n_testing) %dopar% GW_LRC_nonstationary_boot2(dat=dat[,iwk:(n_training_ini - 1 + iwk)]
                                                                         ,fixed_comp=fixed_comp[,iwk:(n_training_ini - 1+ iwk)], 
                                                                         alpha_transf=alpha_transf,fh=(fh-iwk+1), B=B, level=level, 
                                                                         fmethod=fmethod,ncomp_selection =ncomp_selection)
   
    lb <- build.tensor(l=boot,nrows=nrow(sim_data),ncols=fh,name="lb",fh)
    ub <- build.tensor(l=boot,nrows=nrow(sim_data),ncols=fh,name="ub",fh)

     

  }

  result <-compare.data(lb=lb, ub=ub, real=test_data, fh=fh, level)
  
  return(list(lb=lb,ub=ub,result=result))
  
}