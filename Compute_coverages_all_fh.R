library(doMC)
source("CoDa_nonparametric_boot.R")
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



coverage_CoDA <- function(dat,fixed_comp,alpha_transf,fh, B, level, fmethod ="arima",sample_number=n_year,no_core)
{
  sim_data=dat
  grid_point = seq(0, 1, length.out = nrow(sim_data))
  colnames(sim_data) = 1:sample_number
  rownames(sim_data) = grid_point
  
  
  
  # define testing and training samples
  n_val = ncol(sim_data)
  n_testing = round(n_val/5)
  n_training_ini = round(n_val*4/5)
  
  # Construct the testing data
  d_x_t_star_recon_1=GW_LRC_nonstationary_boot(dat,fixed_comp, alpha_transf,fh=fh, B, level, fmethod)$data_untransformed
  test_data = as.matrix(d_x_t_star_recon_1[,(n_training_ini + 1):n_val])
  
  if((n_training_ini + n_testing) != n_val)
  {
    warning("length of training sample + testing sample != total sample")
  }else {
    
    # Computethe lb and up at each forecast horizon.
    
    pointwise_coverage=pointwise_score=CPD_pointwise=matrix(NA,nrow=fh,ncol = 1)

    for (ij in 1:fh) {
      
    pointwise_PI_lb = pointwise_PI_ub = matrix(NA, nrow(sim_data), n_testing)
    registerDoMC(no_core)
    boot = foreach(iwk = 1:n_testing) %dopar% GW_LRC_nonstationary_boot(dat=dat[,iwk:(n_training_ini - ij + iwk)]
                                                                        ,fixed_comp=fixed_comp[,iwk:(n_training_ini - ij + iwk)], 
                                                                        alpha_transf=alpha_transf,fh=ij, B=B, level=level, fmethod=fmethod)
    for(ik in 1:n_testing)
    {
      # pointwise
      pointwise_PI_lb[,ik] = boot[[ik]]$lb[,1]
      pointwise_PI_ub[,ik] = boot[[ik]]$ub[,1]
      
    }
    # compute pointwise coverage probability
    lb_ind = which(pointwise_PI_lb >= test_data)
    ub_ind = which(pointwise_PI_ub <= test_data)
    pointwise_coverage[ij,] = 1 - length(union(lb_ind, ub_ind))/length(test_data)
    
    
    # compute pointwise interval scores
    
    pointwise_score[ij,] = mean(interval_score(holdout = test_data, lb = pointwise_PI_lb,
                                               ub = pointwise_PI_ub, alpha = 0.2),na.rm=TRUE)
    
    CPD_pointwise[ij,]=abs(pointwise_coverage[ij,]-level/100) 
    print(ij)
    }
  
  
  }
  
  result = as.matrix(cbind(pointwise_coverage,pointwise_score,CPD_pointwise))
  colnames(result) = c(paste("Pointwise",level,"coverage"),paste(level,"%", "interval score"),paste("CPD Pointwise", level, "% coverage"))
  
  
  return(result)
  
}