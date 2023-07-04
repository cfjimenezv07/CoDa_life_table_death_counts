# Case: USA

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

##########################
# set a working directory
##########################
setwd("~/My Drive/Summer 2023/STAT 397/PhD project 5/Rcodes/")

##########################
# set a images directory
##########################
dir.l<-"~/My Drive/Summer 2023/STAT 397/PhD project 5/Plots_paper/"
################################################################################
# Dataset entries
################################################################################
year = 1959:2020
n_year = length(year)
age = 0:110
n_age = length(age)
n_prefectures=51

# Row partition
part_list = list()
for(ik in 1:n_prefectures) {
  part_list[[ik]] = (n_year*ik-(n_year-1)):(n_year*ik)
}

#Column partition
n_populations=2
part_list_c = list()
for(ik in 1:n_populations) {
  part_list_c[[ik]] = (n_age*ik-(n_age-1)):(n_age*ik)
}

savepdf <- function(file, width=16, height=10)
{
  fname <- paste(dir.l,file,".pdf",sep="")
  pdf(fname, width=width/2.54, height=height/2.54,
      pointsize=10)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}

################################################################################
# Dataset 
################################################################################
All_USA_male_qx   <- readRDS("./All_USA_male_qx.rds")
All_USA_female_qx <- readRDS("./All_USA_female_qx.rds")

################################################################################
# Get the averages across years for each of the 95 departments the rainbow plot
################################################################################
all_male<-t(list.cbind(All_USA_male_qx))
all_female<-t(list.cbind(All_USA_female_qx))
Male <- matrix(0,nrow = n_age, ncol = n_year)
Female <- matrix(0,nrow = n_age, ncol = n_year)
for (j in 1:n_year) {
  Male[,j]<-colMeans(all_male[seq(from=j,to=n_year*n_prefectures,by=n_year),]) 
  Female[,j]<-colMeans(all_female[seq(from=j,to=n_year*n_prefectures,by=n_year),])
}

################
# rainbow plots
################
par(mfrow=c(1,1))
colnames(Male)<-1:dim(Male)[2]
curves_male<-rainbow::fts(x=1:n_age, y=Male , xname='Age',yname='Life-table death counts')
plot.fds(curves_male,main="USA: male data (1959-2020)")
colnames(Female)<-1:dim(Female)[2]
curves_female<-rainbow::fts(x=1:n_age, y=Female , xname='Age',yname='Life-table death counts')
plot.fds(curves_female,main="USA: female data (1959-2020)")

################################################
#Plots for Figure 1a, 1b.
################################################

#for a specific state
savepdf("Rhode_Island_female")
Female<-All_USA_female_qx[[40]]
colnames(Female)<-1:dim(Female)[2]
curves_female<-rainbow::fts(x=1:n_age, y=Female,xname = "" ,yname='Life-table death counts')
plot.fds(curves_female,ylim=c(0,5200),main="Female")
dev.off()



savepdf("Rhode_Island_male")
Male <- All_USA_male_qx[[40]]
colnames(Male)<-1:dim(Male)[2]
curves_male<-rainbow::fts(x=1:n_age, y=Male,xname = '',yname = '')
plot.fds(curves_male,ylim=c(0,5200),main="Male")
dev.off()


# Compute the Gini-Coefficient. 
require(reldist)
Gini_female<- matrix(0,ncol=n_prefectures,nrow=n_year)
Gini_male<- matrix(0,ncol=n_prefectures,nrow=n_year)
for (i in 1: n_prefectures) {
  Gini_female[,i]=ts(apply(All_USA_female_qx[[i]], 2, gini), start = 1959, end = 2020)
  Gini_male[,i]=ts(apply(All_USA_male_qx[[i]], 2, gini), start = 1959, end = 2020)
}




#####################################################
# Apply the log-ratio transformation
#####################################################
source("CoDa_transformations.R")

Unconstrained_male   <- list()
Unconstrained_female <- list()
alpha_male           <- list()
alpha_female         <- list()
for (i in 1:n_prefectures) {
  transformation <- log_ratio_trans(dat_1=t(All_USA_male_qx[[i]])+1e-5,dat_2=t(All_USA_female_qx[[i]])+1e-5)
  Unconstrained_male[[i]]   <- t(transformation$h_x_t_1)
  Unconstrained_female[[i]] <- t(transformation$h_x_t_2)
  alpha_male[[i]]           <- (transformation$alpha_x_1)
  alpha_female[[i]]         <- (transformation$alpha_x_2)
}

all_unconstrained_male<-t(list.cbind(Unconstrained_male))
all_unconstrained_female<-t(list.cbind(Unconstrained_female))
all_alpha_male   <- list.cbind(alpha_male)
all_alpha_female <- list.cbind(alpha_female)

#####################################################
# FMP decomposition to the unconstrained FTS
#####################################################
source("FMP-ANOVA_decomposition.R")
Y=cbind(all_unconstrained_male,all_unconstrained_female)
Both<-Two_way_median_polish(Y,
                            row_partition_index=part_list,
                            column_partition_index=part_list_c)
Residuals<-Two_way_Residuals(Y,Both,n_prefectures,n_year,n_age)
Res1=Residuals$residuals1
Res2=Residuals$residuals2
Residuals_<-cbind(Res1,Res2)
# Reconstructed data
RR<-Residuals$rd #Matrix with the original data reconstructed from the FMP decomposition
#  It's the proof of the reconstruction of the residuals. 
Residuals$R #The result should be a vector with two entries TRUE, TRUE.
#Indicating that after adding both deterministic and time-varying components the FTS are recovered.
Fixed_part<-Residuals$Fixed_comp # deterministic components to be added up after forecasting

################################################################################
# Computation of the point forecasts
################################################################################
source("Point_forecasting.R")

# With a rolling window approach
cl <- makePSOCKcluster(detectCores()-2)
registerDoParallel(cl)
forecasted_curves_triangular_rolling <- foreach(i = 1:n_prefectures, .packages = c("ftsa")) %dopar% ForecastC(i,fixed_com=Fixed_part,Residuals_f=Residuals_)
stopCluster(cl)


# # With an expanding window approach
cl <- makePSOCKcluster(detectCores()-2)
registerDoParallel(cl)
forecasted_curves_triangular_expanding <- foreach(i = 1:n_prefectures, .packages = c("ftsa")) %dopar% ForecastC_expanding(i,fixed_com=Fixed_part,Residuals_f=Residuals_)
stopCluster(cl)



# Based on the rolling window approach
forecasted_curves        <- list.cbind(t(forecasted_curves_triangular_rolling))
forecasted_curves_male   <- forecasted_curves[1:n_age,]
forecasted_curves_female <- forecasted_curves[(n_age+1):(2*n_age),]


################################################################################
# Transform back to the constrained space
################################################################################
fh=10
dat_transformed_male   <- list()
dat_transformed_female <- list()
Forc_transformed_male   <- list()
Forc_transformed_female <- list()
for (i in 1:n_prefectures) {
  # all the data for each prefecture
  
  alpha_x_1  <- all_alpha_male[,i]
  alpha_x_2  <- all_alpha_female[,i]
  dat_1      <- t(All_USA_male_qx[[i]])+1e-5
  dat_2      <- t(All_USA_female_qx[[i]])+1e-5
  recon_1    <- t(Y[(1+n_year*(i-1)):(n_year*i),1:n_age])
  recon_2    <- t(Y[(1+n_year*(i-1)):(n_year*i),(n_age+1):(2*n_age)])
  fore_val_1 <- forecasted_curves_male[,(fh*i-(fh-1)):(fh*i)]
  fore_val_2 <- forecasted_curves_female[,(fh*i-(fh-1)):(fh*i)]
  
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
  
  dat_transformed_male[[i]]   <- d_x_t_star_recon_1
  dat_transformed_female[[i]] <- d_x_t_star_recon_2

  # back-transformation
  
  f_x_t_star_fore_1 = f_x_t_star_fore_2 = 
    d_x_t_star_fore_1 = d_x_t_star_fore_2 = matrix(NA, n_age, fh)
  for(ik in 1:fh)
  {
    f_x_t_star_fore_1[,ik] = exp(fore_val_1[,ik])/sum(exp(fore_val_1[,ik]))
    f_x_t_star_fore_2[,ik] = exp(fore_val_2[,ik])/sum(exp(fore_val_2[,ik]))
    d_x_t_star_fore_1[,ik] = (f_x_t_star_fore_1[,ik] * alpha_x_1)/sum((f_x_t_star_fore_1[,ik] * alpha_x_1))
    d_x_t_star_fore_2[,ik] = (f_x_t_star_fore_2[,ik] * alpha_x_2)/sum((f_x_t_star_fore_2[,ik] * alpha_x_2))
  }
  colnames(d_x_t_star_fore_1) = colnames(d_x_t_star_fore_2) = 1:fh
  rownames(d_x_t_star_fore_1) = rownames(d_x_t_star_fore_2) = age
  # 
  Forc_transformed_male[[i]] <- d_x_t_star_fore_1
  Forc_transformed_female[[i]] <- d_x_t_star_fore_2
}


####################################################################################
#Compute the error for the point forecasts obtained in the rolling window approach.
####################################################################################
source("forecast_errors.R")
max_h=10
error_MAPE_male<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_KLD_male<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_MAPE_female<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_KLD_female<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_JS1_male<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_JS2_male<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_JS1_female<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_JS2_female<-matrix(0,nrow=max_h,ncol = n_prefectures)
for (i in 1:n_prefectures) {
  male <- dat_transformed_male[[i]]
  female <- dat_transformed_female[[i]]
  forcasted_male<-(Forc_transformed_male[[i]])
  forcasted_female<- (Forc_transformed_female[[i]])
  True_pop_male=(male[,(n_year-(max_h-1)):(n_year)])
  True_pop_female=(female[,(n_year-(max_h-1)):(n_year)])
  for (j in 1:max_h) {
    error_MAPE_male[j,i]<-mape(forecast=forcasted_male[,j], true=True_pop_male[,j])
    error_KLD_male[j,i]<-KLD(forecast=forcasted_male[,j], true=True_pop_male[,j])
    error_MAPE_female[j,i]<-mape(forecast=forcasted_female[,j], true=True_pop_female[,j])
    error_KLD_female[j,i]<-KLD(forecast=forcasted_female[,j], true=True_pop_female[,j])
    error_JS1_male[j,i]<-JSD_simple(forecast=forcasted_male[,j], true=True_pop_male[,j])
    error_JS2_male[j,i]<-JSD_geom(forecast=forcasted_male[,j], true=True_pop_male[,j])
    error_JS1_female[j,i]<-JSD_simple(forecast=forcasted_female[,j], true=True_pop_female[,j])
    error_JS2_female[j,i]<-JSD_geom(forecast=forcasted_female[,j], true=True_pop_female[,j])
  }
  
}

All_errors_rolling<-list(error_MAPE_male,error_MAPE_female,
                         error_KLD_male,error_KLD_female,
                         error_JS1_male,error_JS1_female,
                         error_JS2_male,error_JS2_female)

Errors_mean_rolling<-matrix(0,nrow=n_prefectures,ncol=length(All_errors_rolling))
Errors_sd_rolling<-matrix(0,nrow=n_prefectures,ncol=length(All_errors_rolling))
for (i in 1:length(All_errors_rolling)) {
  error=All_errors_rolling[[i]]
  Errors_mean_rolling[,i]=apply( error,2,mean)
  Errors_sd_rolling[,i]=apply( error,2,sd)
}
colnames(Errors_mean_rolling)<-c("MAPE_male","MAPE_female","KLD_male","KLD_female","JS1_male","JS1_female","JS2_male","JS2_female")
colnames(Errors_sd_rolling)<-c("MAPE_male","MAPE_female","KLD_male","KLD_female","JS1_male","JS1_female","JS2_male","JS2_female")
rownames(Errors_mean_rolling)<-names_states
rownames(Errors_sd_rolling)<-names_states


########################################################################
# Compute interval forecasts and coverage probability
########################################################################
# Split residuals per prefecture
male_prefecture_res<-lapply(1:length(part_list), 
                            function(k){Res1[part_list[[k]], ]})
female_prefecture_res<-lapply(1:length(part_list), 
                              function(k){Res2[part_list[[k]], ]})

#split the fixed components by prefecture
male_prefecture_fixed<-lapply(1:length(part_list), 
                              function(k){Fixed_part[part_list[[k]],1:n_age ]})

female_prefecture_fixed<-lapply(1:length(part_list), 
                                function(k){Fixed_part[part_list[[k]],(n_age+1):(2*n_age) ]})

source("Compute_coverages_boot.R")
Emp_cov_male_80<-matrix(0,nrow=n_prefectures,ncol = 3)
Emp_cov_male_95<-matrix(0,nrow=n_prefectures,ncol = 3)
Emp_cov_female_80<-matrix(0,nrow=n_prefectures,ncol = 3)
Emp_cov_female_95<-matrix(0,nrow=n_prefectures,ncol = 3)
for (i in 1:n_prefectures) {
  Emp_cov_male_80[i,]<-coverage_CoDA(dat=t(male_prefecture_res[[i]]),fixed_comp=t(male_prefecture_fixed[[i]]),alpha_transf=all_alpha_male[,i],
                                     sample_number=n_year, fh=10, B=1000, level=80, fmethod ="arima",no_core=detectCores()-2)
  Emp_cov_male_95[i,]<-coverage_CoDA(dat=t(male_prefecture_res[[i]]),fixed_comp=t(male_prefecture_fixed[[i]]),alpha_transf=all_alpha_male[,i],
                                     sample_number=n_year, fh=10, B=1000, level=95, fmethod ="arima",no_core=detectCores()-2)
  Emp_cov_female_80[i,]<-coverage_CoDA(dat=t(female_prefecture_res[[i]]),fixed_comp=t(female_prefecture_fixed[[i]]),alpha_transf=all_alpha_female[,i],
                                       sample_number=n_year, fh=10, B=1000, level=80, fmethod ="arima",no_core=detectCores()-2)
  Emp_cov_female_95[i,]<-coverage_CoDA(dat=t(female_prefecture_res[[i]]),fixed_comp=t(female_prefecture_fixed[[i]]),alpha_transf=all_alpha_female[,i],
                                       sample_number=n_year, fh=10, B=1000, level=95, fmethod ="arima",no_core=detectCores()-2)
}

Emp_cov_80<-cbind(Emp_cov_male_80,Emp_cov_female_80)
colnames(Emp_cov_80) = c("Pointwise 80% coverage",
                         "80% interval score",
                         "CPD Pointwise 80% coverage",
                         "Pointwise 80% coverage", 
                         "80% interval score",
                         "CPD Pointwise 80% coverage")
rownames(Emp_cov_80)<-names_states
Emp_cov_95<-cbind(Emp_cov_male_95,Emp_cov_female_95)
colnames(Emp_cov_95) = c("Pointwise 95% coverage",
                         "95% interval score",
                         "CPD Pointwise 95% coverage",
                         "Pointwise 95% coverage", 
                         "95% interval score",
                         "CPD Pointwise 95% coverage")
rownames(Emp_cov_95)<-names_states


########################################################################
# Comparison with the functional mean ANOVA approach (FM-ANOVA)
########################################################################
source("FM-ANOVA_decomposition.R")

# This function computes the functional mean ANOVA decomposition based on means
FANOVA_means<-FANOVA(data_pop1=t(all_unconstrained_male),data_pop2=t(all_unconstrained_female),n_year
                     ,n_prefectures,n_age,n_populations=2,row_par=part_list)
#This function computes the functional residuals after removing the deterministic components
# obtained from the FANOVA function.
Residuals_means<-Two_way_Residuals_means(FANOVA_means,data_pop1=t(all_unconstrained_male),data_pop2=t(all_unconstrained_female)
                                         ,n_prefectures,n_year,n_age)
Res1_means=Residuals_means$residuals1_mean
Res2_means=Residuals_means$residuals2_mean
Residuals_mean<-cbind(Res1_means,Res2_means)
# Reconstructed data
RR<-Residuals_means$rd #Matrix with the original data reconstructed from the FMP decomposition
#  It's the proof of the reconstruction of the residuals. 
Residuals_means$R #The result should be a vector with two entries TRUE, TRUE.
#Indicating that after adding both deterministic and time-varying components the FTS are recovered.
Fixed_part_means<-Residuals_means$Fixed_comp_mean # deterministic components to be added up after forecasting


################################################################################
# Computation of the point forecasts based on functional mean ANOVA (FM-ANOVA)
################################################################################
source("Point_forecasting_FANOVA.R")
# With a rolling window approach
cl <- makePSOCKcluster(detectCores()-2)
registerDoParallel(cl)
forecasted_curves_triangular_means_rolling <- foreach(i = 1:n_prefectures, .packages = c("ftsa")) %dopar% ForecastC_mean(i,fixed_com=Fixed_part_means,Residuals_f=Residuals_mean)
stopCluster(cl)

 
# # With an expanding window approach
cl <- makePSOCKcluster(detectCores()-2)
registerDoParallel(cl)
forecasted_curves_triangular_means_expanding <- foreach(i = 1:n_prefectures, .packages = c("ftsa")) %dopar% ForecastC_mean_expanding(i,fixed_com=Fixed_part_means,Residuals_f=Residuals_mean)
stopCluster(cl)


# Based on rolling window
forecasted_curves_triangular_means_rolling <- readRDS("./forecasted_curves_triangular_means_rolling.rds")
forecasted_curves_mean        <- list.cbind(t(forecasted_curves_triangular_means_rolling))
forecasted_curves_mean_male   <- forecasted_curves_mean[1:n_age,]
forecasted_curves_mean_female <- forecasted_curves_mean[(n_age+1):(2*n_age),]


################################################################################
# Transform back to the constrained space
################################################################################
fh=10
dat_transformed_mean_male   <- list()
dat_transformed_mean_female <- list()
Forc_transformed_mean_male   <- list()
Forc_transformed_mean_female <- list()
for (i in 1:n_prefectures) {
  # all the data for each prefecture
  
  alpha_x_1  <- all_alpha_male[,i]
  alpha_x_2  <- all_alpha_female[,i]
  dat_1      <- t(All_USA_male_qx[[i]])+1e-5
  dat_2      <- t(All_USA_female_qx[[i]])+1e-5
  recon_1    <- t(Y[(1+n_year*(i-1)):(n_year*i),1:n_age])
  recon_2    <- t(Y[(1+n_year*(i-1)):(n_year*i),(n_age+1):(2*n_age)])
  fore_val_1 <- forecasted_curves_mean_male[,(fh*i-(fh-1)):(fh*i)]
  fore_val_2 <- forecasted_curves_mean_female[,(fh*i-(fh-1)):(fh*i)]
  
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
  
  dat_transformed_mean_male[[i]]   <- d_x_t_star_recon_1
  dat_transformed_mean_female[[i]] <- d_x_t_star_recon_2
 
  # back-transformation
  
  f_x_t_star_fore_1 = f_x_t_star_fore_2 = 
    d_x_t_star_fore_1 = d_x_t_star_fore_2 = matrix(NA, n_age, fh)
  for(ik in 1:fh)
  {
    f_x_t_star_fore_1[,ik] = exp(fore_val_1[,ik])/sum(exp(fore_val_1[,ik]))
    f_x_t_star_fore_2[,ik] = exp(fore_val_2[,ik])/sum(exp(fore_val_2[,ik]))
    d_x_t_star_fore_1[,ik] = (f_x_t_star_fore_1[,ik] * alpha_x_1)/sum((f_x_t_star_fore_1[,ik] * alpha_x_1))
    d_x_t_star_fore_2[,ik] = (f_x_t_star_fore_2[,ik] * alpha_x_2)/sum((f_x_t_star_fore_2[,ik] * alpha_x_2))
  }
  colnames(d_x_t_star_fore_1) = colnames(d_x_t_star_fore_2) = 1:fh
  rownames(d_x_t_star_fore_1) = rownames(d_x_t_star_fore_2) = age
  # 
  Forc_transformed_mean_male[[i]] <- d_x_t_star_fore_1
  Forc_transformed_mean_female[[i]] <- d_x_t_star_fore_2
}

##############################################################################################################
#Compute the error based on FM-ANOVA approach, with the forecasts obtained in the rolling window approach
##############################################################################################################
source("forecast_errors.R")
max_h=10
error_MAPE_means_male<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_KLD_means_male<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_MAPE_means_female<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_KLD_means_female<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_JS1_means_male<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_JS2_means_male<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_JS1_means_female<-matrix(0,nrow=max_h,ncol = n_prefectures)
error_JS2_means_female<-matrix(0,nrow=max_h,ncol = n_prefectures)
for (i in 1:n_prefectures) {
  male <- dat_transformed_mean_male [[i]]
  female <- dat_transformed_mean_female [[i]]
  forecasted_male<-(Forc_transformed_mean_male[[i]])
  forecasted_female<- (Forc_transformed_mean_female[[i]])
  True_pop_male=(male[,(n_year-(max_h-1)):(n_year)])
  True_pop_female=(female[,(n_year-(max_h-1)):(n_year)])
  for (j in 1:max_h) {
    error_MAPE_means_male[j,i]<-mape(forecast=forecasted_male[,j], true=True_pop_male[,j])
    error_KLD_means_male[j,i]<-KLD(forecast=forecasted_male[,j], true=True_pop_male[,j])
    error_MAPE_means_female[j,i]<-mape(forecast=forecasted_female[,j], true=True_pop_female[,j])
    error_KLD_means_female[j,i]<-KLD(forecast=forecasted_female[,j], true=True_pop_female[,j])
    error_JS1_means_male[j,i]<-JSD_simple(forecast=forecasted_male[,j], true=True_pop_male[,j])
    error_JS2_means_male[j,i]<-JSD_geom(forecast=forecasted_male[,j], true=True_pop_male[,j])
    error_JS1_means_female[j,i]<-JSD_simple(forecast=forecasted_female[,j], true=True_pop_female[,j])
    error_JS2_means_female[j,i]<-JSD_geom(forecast=forecasted_female[,j], true=True_pop_female[,j])
  }
  
}

All_errors_mean_rolling<-list(error_MAPE_means_male,error_MAPE_means_female,
                              error_KLD_means_male,error_KLD_means_female,
                              error_JS1_means_male,error_JS1_means_female,
                              error_JS2_means_male,error_JS2_means_female)
Errors_mean_basedmeans_rolling<-matrix(0,nrow=n_prefectures,ncol=length(All_errors_mean_rolling))
Errors_sd_basedmeans_rolling<-matrix(0,nrow=n_prefectures,ncol=length(All_errors_mean_rolling))
for (i in 1:length(All_errors_mean_rolling)) {
  error=All_errors_mean_rolling[[i]]
  Errors_mean_basedmeans_rolling[,i]=apply(error,2,mean)
  Errors_sd_basedmeans_rolling[,i]=apply( error,2,sd)
}
colnames(Errors_mean_basedmeans_rolling)<-c("MAPE_male","MAPE_female","KLD_male","KLD_female","JS1_male","JS1_female","JS2_male","JS2_female")
colnames(Errors_sd_basedmeans_rolling)<-c("MAPE_male","MAPE_female","KLD_male","KLD_female","JS1_male","JS1_female","JS2_male","JS2_female")
rownames(Errors_mean_basedmeans_rolling)<-names_states
rownames(Errors_sd_basedmeans_rolling)<-names_states

############################################################################
# Compute interval forecasts and coverage probability for FANOVA approach
############################################################################
# Split residuals per prefecture
male_prefecture_res_means<-lapply(1:length(part_list), 
                                  function(k){Res1_means[part_list[[k]], ]})
female_prefecture_res_means<-lapply(1:length(part_list), 
                                    function(k){Res2_means[part_list[[k]], ]})

#split the fixed components by prefecture
male_prefecture_fixed_means<-lapply(1:length(part_list), 
                                    function(k){Fixed_part_means[part_list[[k]],1:n_age ]})

female_prefecture_fixed_means<-lapply(1:length(part_list), 
                                      function(k){Fixed_part_means[part_list[[k]],(n_age+1):(2*n_age) ]})

source("Compute_coverages_boot.R")
Emp_cov_male_80_means<-matrix(0,nrow=n_prefectures,ncol = 3)
Emp_cov_male_95_means<-matrix(0,nrow=n_prefectures,ncol = 3)
Emp_cov_female_80_means<-matrix(0,nrow=n_prefectures,ncol = 3)
Emp_cov_female_95_means<-matrix(0,nrow=n_prefectures,ncol = 3)
for (i in 1:n_prefectures) {
  Emp_cov_male_80_means[i,]<-coverage_CoDA(dat=t(male_prefecture_res_means[[i]]),fixed_comp=t(male_prefecture_fixed_means[[i]]),alpha_transf=all_alpha_male[,i],
                                           sample_number=n_year, fh=10, B=1000, level=80, fmethod ="arima",no_core=detectCores()-2)
  Emp_cov_male_95_means[i,]<-coverage_CoDA(dat=t(male_prefecture_res_means[[i]]),fixed_comp=t(male_prefecture_fixed_means[[i]]),alpha_transf=all_alpha_male[,i],
                                           sample_number=n_year, fh=10, B=1000, level=95, fmethod ="arima",no_core=detectCores()-2)
  Emp_cov_female_80_means[i,]<-coverage_CoDA(dat=t(female_prefecture_res_means[[i]]),fixed_comp=t(female_prefecture_fixed_means[[i]]),alpha_transf=all_alpha_female[,i],
                                             sample_number=n_year, fh=10, B=1000, level=80, fmethod ="arima",no_core=detectCores()-2)
  Emp_cov_female_95_means[i,]<-coverage_CoDA(dat=t(female_prefecture_res_means[[i]]),fixed_comp=t(female_prefecture_fixed_means[[i]]),alpha_transf=all_alpha_female[,i],
                                             sample_number=n_year, fh=10, B=1000, level=95, fmethod ="arima",no_core=detectCores()-2)

}

Emp_cov_80_means<-cbind(Emp_cov_male_80_means,Emp_cov_female_80_means)
colnames(Emp_cov_80_means) = c("Pointwise 80% coverage",
                               "80% interval score",
                               "CPD Pointwise 80% coverage",
                               "Pointwise 80% coverage", 
                               "80% interval score",
                               "CPD Pointwise 80% coverage")
rownames(Emp_cov_80_means)<-names_states
Emp_cov_95_means<-cbind(Emp_cov_male_95_means,Emp_cov_female_95_means)
colnames(Emp_cov_95_means) = c("Pointwise 95% coverage",
                               "95% interval score",
                               "CPD Pointwise 95% coverage",
                               "Pointwise 95% coverage", 
                               "95% interval score",
                               "CPD Pointwise 95% coverage")
rownames(Emp_cov_95_means)<-names_states

