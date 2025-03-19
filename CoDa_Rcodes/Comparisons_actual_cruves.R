#Comparison of forecast curves with TNH, GSY, and PLS methods. 
setwd("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/Github/CoDa_life_table_death_counts/CoDa_Rcodes")
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/Github/CoDa_life_table_death_counts/CoDa_Rcodes/Results/"
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


################################################################################
# Dataset entries
################################################################################
dataset_entries<-readRDS("./dataset_entries/USA/dataset_entries.rds")
year = dataset_entries[[1]]
n_year = dataset_entries[[2]]
age = dataset_entries[[3]]
n_age = dataset_entries[[4]]
n_states=dataset_entries[[5]]



# Row partition
part_list = list()
for(ik in 1:n_states) {
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
names_states       <- readRDS("./datasets/USA/names_states.rds")
All_USA_male_qx    <- readRDS("./datasets/USA/All_USA_male_qx.rds")
All_USA_female_qx  <- readRDS("./datasets/USA/All_USA_female_qx.rds")

################################################################################
# Get the averages across years for each of the 95 departments the rainbow plot
################################################################################
all_male<-t(list.cbind(All_USA_male_qx))
all_female<-t(list.cbind(All_USA_female_qx))
Male <- matrix(0,nrow = n_age, ncol = n_year)
Female <- matrix(0,nrow = n_age, ncol = n_year)
for (j in 1:n_year) {
  Male[,j]<-colMeans(all_male[seq(from=j,to=n_year*n_states,by=n_year),]) 
  Female[,j]<-colMeans(all_female[seq(from=j,to=n_year*n_states,by=n_year),])
}

##################################################################
# 1.  FMP-ANOVA + FM-ANOVA
##################################################################
#####################################################
# Apply the log-ratio transformation
#####################################################
source("CoDa_transformations.R")

Unconstrained_male   <- list()
Unconstrained_female <- list()
alpha_male           <- list()
alpha_female         <- list()
for (i in 1:n_states) {
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
library("ftsa")
Y=cbind(all_unconstrained_male,all_unconstrained_female)
Both<-ftsa::Two_way_median_polish(Y,year,age,n_states,n_populations)
Residuals<- ftsa::Two_way_Residuals(Y,n_states,year,age,n_populations)
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
source("New_Point_forecast2.R")

# With a rolling window approach
cl <- makePSOCKcluster(detectCores()-2)
registerDoParallel(cl)
forecasted_FMP_lrc_ARIMA_USA_3 <- foreach(i = 1:n_states, .packages = c("ftsa","doMC")) %do% {ForecastC(i,fixed_com=Fixed_part,Residuals_f=Residuals_,est_method = "lrc",
                                                                                                        prediction_method="ARIMA",select_K="EVR", K=6,no_core=detectCores()-2,
                                                                                                        all_alpha_male= all_alpha_male,all_alpha_female=all_alpha_female,
                                                                                                        n_year=n_year,n_age=n_age,age=age)}

forecasted_FMP_lrc_ARIMA_USA_4 <- foreach(i = 1:n_states, .packages = c("ftsa","doMC")) %do% {ForecastC(i,fixed_com=Fixed_part,Residuals_f=Residuals_,est_method = "lrc",
                                                                                                        prediction_method="ARIMA",select_K="Fixed",K=6,no_core=detectCores()-2,
                                                                                                        all_alpha_male= all_alpha_male,all_alpha_female=all_alpha_female,
                                                                                                        n_year=n_year,n_age=n_age,age=age)}
stopCluster(cl)




# saveRDS(forecasted_FMP_lrc_ARIMA_USA_3,paste0(dirl.p,"forecasted_FMP_lrc_ARIMA_USA_3.rds"))
# saveRDS(forecasted_FMP_lrc_ARIMA_USA_4,paste0(dirl.p,"forecasted_FMP_lrc_ARIMA_USA_4.rds"))


####################################################################################
#Transform the forecasted curves_back
####################################################################################
# EVR
dat_transformed_female_EVR_FMP  <- list()
dat_transformed_male_EVR_FMP    <- list()
Forc_transformed_female_EVR_FMP <- list()
Forc_transformed_male_EVR_FMP   <- list()
for (i in 1:n_states) {
  
  # forecasted _3 if EVR and _4 if K=6
  forecasted <-forecasted_FMP_lrc_ARIMA_USA_3
  
  dat_transformed_male_EVR_FMP[[i]] <- forecasted[[i]]$dat_transformed_male
  dat_transformed_female_EVR_FMP[[i]] <- forecasted[[i]]$dat_transformed_female
  
  # transform forecast
  Forc_transformed_male_EVR_FMP[[i]]   <- forecasted[[i]]$Forc_transformed_male
  Forc_transformed_female_EVR_FMP[[i]] <- forecasted[[i]]$Forc_transformed_female
  
}

saveRDS(dat_transformed_female_EVR_FMP,paste0(dir.r,"dat_transformed_female_EVR_FMP.rds"))
saveRDS(dat_transformed_male_EVR_FMP,paste0(dir.r,"dat_transformed_male_EVR_FMP.rds"))
saveRDS(Forc_transformed_female_EVR_FMP,paste0(dir.r,"Forc_transformed_female_EVR_FMP.rds"))
saveRDS(Forc_transformed_male_EVR_FMP,paste0(dir.r,"Forc_transformed_male_EVR_FMP.rds"))

# K=6
dat_transformed_female_K_FMP  <- list()
dat_transformed_male_K_FMP    <- list()
Forc_transformed_female_K_FMP <- list()
Forc_transformed_male_K_FMP   <- list()
for (i in 1:n_states) {
  
  # forecasted _3 if K and _4 if K=6
  forecasted <-forecasted_FMP_lrc_ARIMA_USA_4
  
  dat_transformed_male_K_FMP[[i]] <- forecasted[[i]]$dat_transformed_male
  dat_transformed_female_K_FMP[[i]] <- forecasted[[i]]$dat_transformed_female
  
  # transform forecast
  Forc_transformed_male_K_FMP[[i]]   <- forecasted[[i]]$Forc_transformed_male
  Forc_transformed_female_K_FMP[[i]] <- forecasted[[i]]$Forc_transformed_female
  
}

saveRDS(dat_transformed_female_K_FMP,paste0(dir.r,"dat_transformed_female_K_FMP.rds"))
saveRDS(dat_transformed_male_K_FMP,paste0(dir.r,"dat_transformed_male_K_FMP.rds"))
saveRDS(Forc_transformed_female_K_FMP,paste0(dir.r,"Forc_transformed_female_K_FMP.rds"))
saveRDS(Forc_transformed_male_K_FMP,paste0(dir.r,"Forc_transformed_male_K_FMP.rds"))


########################################################################
# Comparison with the functional mean ANOVA approach (FM-ANOVA)
########################################################################

# This function computes the functional mean ANOVA decomposition based on means
FANOVA_means<-ftsa::FANOVA(data_pop1=t(all_unconstrained_male),data_pop2=t(all_unconstrained_female),year,age,n_states,n_populations)

#This function computes the functional residuals after removing the deterministic components
# obtained from the FANOVA function.
Residuals_means<-ftsa::Two_way_Residuals_means(data_pop1=t(all_unconstrained_male),data_pop2=t(all_unconstrained_female)
                                               ,year,age,n_states,n_populations)

Res1_means=Residuals_means$residuals1_mean
Res2_means=Residuals_means$residuals2_mean
Residuals_mean<-cbind(Res1_means,Res2_means)

# Reconstructed data
RR<-Residuals_means$rd #Matrix with the original data reconstructed from the FMP decomposition
#  It's the proof of the reconstruction of the residuals. 
Residuals_means$R #The result should be a vector with two entries TRUE, TRUE.
#Indicating that after adding both deterministic and time-varying components the FTS are recovered.
Fixed_part_means<-Residuals_means$Fixed_comp_mean # deterministic components to be added up after forecasting
Fixed_part_means_1 <- Fixed_part_means[,1:n_age]
Fixed_part_means_2 <- Fixed_part_means[,(n_age+1):(2*n_age)]

################################################################################
# Computation of the point forecasts based on functional mean ANOVA (FM-ANOVA)
################################################################################
source("New_Point_forecast2.R")

# With a rolling window approach
cl <- makePSOCKcluster(detectCores()-2)
registerDoParallel(cl)
forecasted_FM_lrc_ARIMA_USA_3 <- foreach(i = 1:n_states, .packages = c("ftsa","doMC")) %do% {ForecastC(i,fixed_com=Fixed_part_means,Residuals_f=Residuals_mean,est_method = "lrc",prediction_method = "ARIMA",select_K="EVR",K=6
                                                                                                       ,no_core=detectCores()-2,
                                                                                                       all_alpha_male= all_alpha_male,all_alpha_female=all_alpha_female,
                                                                                                       n_year=n_year,n_age=n_age,age=age)}
forecasted_FM_lrc_ARIMA_USA_4 <- foreach(i = 1:n_states, .packages = c("ftsa","doMC")) %do% {ForecastC(i,fixed_com=Fixed_part_means,Residuals_f=Residuals_mean,est_method = "lrc",prediction_method = "ARIMA",select_K="Fixed",K=6,no_core=detectCores()-2,
                                                                                                       all_alpha_male= all_alpha_male,all_alpha_female=all_alpha_female,
                                                                                                       n_year=n_year,n_age=n_age,age=age)}
stopCluster(cl)


# saveRDS(forecasted_FM_lrc_ARIMA_USA_3,paste0(dirl.p,"forecasted_FM_lrc_ARIMA_USA_3.rds"))
# saveRDS(forecasted_FM_lrc_ARIMA_USA_4,paste0(dirl.p,"forecasted_FM_lrc_ARIMA_USA_4.rds"))

# EVR
dat_transformed_female_EVR_FM  <- list()
dat_transformed_male_EVR_FM    <- list()
Forc_transformed_female_EVR_FM <- list()
Forc_transformed_male_EVR_FM   <- list()
for (i in 1:n_states) {
  
  # forecasted _3 if EVR and _4 if K=6
  forecasted <-forecasted_FM_lrc_ARIMA_USA_3
  
  dat_transformed_male_EVR_FM[[i]] <- forecasted[[i]]$dat_transformed_male
  dat_transformed_female_EVR_FM[[i]] <- forecasted[[i]]$dat_transformed_female
  
  # transform forecast
  Forc_transformed_male_EVR_FM[[i]]   <- forecasted[[i]]$Forc_transformed_male
  Forc_transformed_female_EVR_FM[[i]] <- forecasted[[i]]$Forc_transformed_female
  
}

saveRDS(dat_transformed_female_EVR_FM,paste0(dir.r,"dat_transformed_female_EVR_FM.rds"))
saveRDS(dat_transformed_male_EVR_FM,paste0(dir.r,"dat_transformed_male_EVR_FM.rds"))
saveRDS(Forc_transformed_female_EVR_FM,paste0(dir.r,"Forc_transformed_female_EVR_FM.rds"))
saveRDS(Forc_transformed_male_EVR_FM,paste0(dir.r,"Forc_transformed_male_EVR_FM.rds"))

# K=6
dat_transformed_female_K_FM  <- list()
dat_transformed_male_K_FM    <- list()
Forc_transformed_female_K_FM <- list()
Forc_transformed_male_K_FM   <- list()
for (i in 1:n_states) {
  
  # forecasted _3 if K and _4 if K=6
  forecasted <-forecasted_FM_lrc_ARIMA_USA_4
  
  dat_transformed_male_K_FM[[i]] <- forecasted[[i]]$dat_transformed_male
  dat_transformed_female_K_FM[[i]] <- forecasted[[i]]$dat_transformed_female
  
  # transform forecast
  Forc_transformed_male_K_FM[[i]]   <- forecasted[[i]]$Forc_transformed_male
  Forc_transformed_female_K_FM[[i]] <- forecasted[[i]]$Forc_transformed_female
  
}

saveRDS(dat_transformed_female_K_FM,paste0(dir.r,"dat_transformed_female_K_FM.rds"))
saveRDS(dat_transformed_male_K_FM,paste0(dir.r,"dat_transformed_male_K_FM.rds"))
saveRDS(Forc_transformed_female_K_FM,paste0(dir.r,"Forc_transformed_female_K_FM.rds"))
saveRDS(Forc_transformed_male_K_FM,paste0(dir.r,"Forc_transformed_male_K_FM.rds"))


################################################################################
# auxiliary functions for the computation of the TNH and GSY methods
################################################################################
source('./MITS_class.R') 
source("./aux_TNH.R")
#####
# GSY
####



