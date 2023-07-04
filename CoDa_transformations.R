 # CoDa transformation algorthim

# Implementation of the centre log-ratio-transformation
##########################
# set a working directory
##########################
setwd("~/My Drive/Summer 2023/STAT 397/PhD project 5/Rcodes/")

require(psych)


# dat_1: first data set (n_year by n_age)  
# dat_2: second data set (n_year by n_age)
# fh: forecast horizon
# n_prefectures: states or departments
# modeling_method: Functional Median Polish (FMP) and Functional Mean ANOVA (FM)
# forecasting_method: random walk without drift, random walk with drift, exponential smoothing, autoregressive integrated moving average

log_ratio_trans <- function(dat_1, dat_2)
{
  n_year = nrow(dat_1)
  n_age = ncol(dat_1)
  year_index = rownames(dat_1)
  age_index = colnames(dat_1)
  
  #1.  standardize life table death to sum to 1
  
  dat_1_center = sweep(dat_1, 1, apply(dat_1, 1, sum), "/")
  dat_2_center = sweep(dat_2, 1, apply(dat_2, 1, sum), "/")
  
  alpha_x_1 = alpha_x_2 = vector("numeric", n_age)
  for(ik in 1:n_age)
  {
    alpha_x_1[ik] = geometric.mean(dat_1_center[,ik])
    alpha_x_2[ik] = geometric.mean(dat_2_center[,ik])
  }
  
  f_x_t_1 = f_x_t_2 = matrix(NA, n_year, n_age)
  for(ik in 1:n_year)
  {
    f_x_t_1[ik,] = (dat_1[ik,]/alpha_x_1)/sum(dat_1[ik,]/alpha_x_1)
    f_x_t_2[ik,] = (dat_2[ik,]/alpha_x_2)/sum(dat_2[ik,]/alpha_x_2)
  }
  
  # 2. Apply the centred log-ratio transformation
  
  g_t_1 = g_t_2 = vector("numeric", n_year)
  h_x_t_1 = h_x_t_2 = matrix(NA, n_year, n_age)
  for(ik in 1:n_year)
  {
    g_t_1[ik] = geometric.mean(f_x_t_1[ik,])
    h_x_t_1[ik,] = log(f_x_t_1[ik,]/g_t_1[ik])
    
    g_t_2[ik] = geometric.mean(f_x_t_2[ik,])
    h_x_t_2[ik,] = log(f_x_t_2[ik,]/g_t_2[ik])
  }
  colnames(h_x_t_1) = colnames(h_x_t_2) = age_index
  rownames(h_x_t_1) = rownames(h_x_t_2) = year_index

    return(list(h_x_t_1 = h_x_t_1, h_x_t_2 = h_x_t_2,alpha_x_1 = alpha_x_1,alpha_x_2 = alpha_x_2))
  
}

