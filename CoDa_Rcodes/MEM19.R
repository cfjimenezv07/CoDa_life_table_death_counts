# Alternative method of PLC19 
#"The maximum entropy mortality model: forecasting mortality using statistical moments"

#"The maximum entropy mortality model: forecasting mortality using statistical moments"
# Replicate Pascariu et. al. 2019 (PLC) MEM method for our datasets. 
# These Rcodes are public posted in https://github.com/mpascariu/MortalityForecast.
# We adapted to our paper case.



###########################################################################################
# 1. First install all the R packages and libraries required for the project
###########################################################################################
packages <- c("generics", "demography", "forecast","fda","fdaoutlier", 
              "rlist", "mrfDepth","ftsa","rainbow")

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


# Mortality forecast package
require(devtools)
devtools::install_github("mpascariu/MortalityForecast")
library("MortalityForecast")
##########################
#2.  set a working directory
##########################
# Set a working directory
setwd("~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@kaust.edu.sa/My Drive/Spring_2024/STAT_397/PhD_Project_5/Codes_paper2/")
##########################
# set a images directory
##########################
dir.l<-"~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@kaust.edu.sa/My Drive/Spring_2024/STAT_397/PhD_Project_5/Plots_paper/"
##########################
# set a results directory
##########################
dirl.p<-"~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@kaust.edu.sa/My Drive/Spring_2024/STAT_397/PhD_Project_5/Resuts_final/USA/"


##########################################################################################################
# 3. Load dataset entries and define the row (by states) and column partitions (by gender)
# Note: If you want to change to another dataset, please change the word USA onwards. (e.g. France/Japan) 
##########################################################################################################
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
################################################################################
# Dataset 
################################################################################
names_states       <- readRDS("./datasets/USA/names_states.rds")
All_USA_male_qx    <- readRDS("./datasets/USA/All_USA_male_qx.rds")
All_USA_female_qx  <- readRDS("./datasets/USA/All_USA_female_qx.rds")


################################################################################
# Forecast from MEM19
################################################################################
# With rolling window approach
################################################################################
max_h = 10
Forecast_female <- list()
Forecast_male   <- list()

for (i in 1:n_states) {
  forecasted_prefecture_female <- array(0, dim=c(n_age,max_h,max_h))
  forecasted_prefecture_male   <- array(0, dim=c(n_age,max_h,max_h))
  data_male   = All_USA_male_qx[[i]]
  data_female = All_USA_female_qx[[i]]
  for (k in 1:max_h) {
    y=year[k:(n_year - max_h - 1 + k)]
    last_year <- y[length(y)]
    data_1 =data_male[,k:(n_year - max_h - 1 + k)]
    data_2 =data_female[,k:(n_year - max_h - 1 + k)]
    male_model_MEM = model.MEM(data =data_1  , x = age , y = y , n = 4)
    female_model_MEM = model.MEM(data = data_2 , x = age , y = y , n = 4)
    year_ranges <- as.character(last_year + 1:(max_h - k +1))
    if (k == 1) {
      forecasted_prefecture_male[,,k] = predict(male_model_MEM, h = max_h, x.h = age)$predicted.values[, year_ranges]*10^5
      forecasted_prefecture_female[,,k] = predict(female_model_MEM, h = max_h, x.h = age)$predicted.values[, year_ranges]*10^5
    } else {
      start <- max_h - k + 2
      forecasted_prefecture_male[,-(start:max_h),k] = predict(male_model_MEM, h = max_h, x.h = age)$predicted.values[, year_ranges]*10^5
      forecasted_prefecture_female[,-(start:max_h),k] = predict(female_model_MEM, h = max_h, x.h = age)$predicted.values[, year_ranges]*10^5
    }

  }
  Forecast_female[[i]] <- forecasted_prefecture_female
  Forecast_male[[i]]   <- forecasted_prefecture_male
  print(i)
}

# saveRDS(Forecast_female,paste0(dirl.p,"Forecast_female_MEM_USA_1.rds"))
# saveRDS(Forecast_male,paste0(dirl.p,"Forecast_male_MEM_USA_1.rds"))
##############################################################################################################
#Compute the point forecast errors obtained in the rolling window approach for  MEM method.
##############################################################################################################
source("./forecast_errors.R")
h_max=10
# TNH
MEM_error_KLD_male<-matrix(0,nrow=h_max,ncol = n_states)
MEM_error_KLD_female<-matrix(0,nrow=h_max,ncol = n_states)
MEM_error_JS1_male<-matrix(0,nrow=h_max,ncol = n_states)
MEM_error_JS1_female<-matrix(0,nrow=h_max,ncol = n_states)
MEM_error_JS2_male<-matrix(0,nrow=h_max,ncol = n_states)
MEM_error_JS2_female<-matrix(0,nrow=h_max,ncol = n_states)


for (i in 1:n_states) {
  
  
  dat_male   = All_USA_male_qx[[i]]
  dat_female = All_USA_female_qx[[i]]
  
  
  forecast_female <-Forecast_female[[i]]
  forecast_male   <-Forecast_male[[i]]
  
  for (j in 1:h_max) {
    
    
    errors_female  <- rowMeans(sapply( 1:(h_max-j+1),pairwise_errors, forecast_female[,j,],dat_female[,(n_year-h_max+j):(n_year)])) 
    errors_male  <- rowMeans(sapply( 1:(h_max-j+1),pairwise_errors,  forecast_male[,j,],dat_male[,(n_year-h_max+j):(n_year)])) 
    
    
    
    MEM_error_KLD_female[j, i]  <-  errors_female[2]
    MEM_error_JS1_female[j, i]  <-  errors_female[3]
    MEM_error_JS2_female[j, i]  <-  errors_female[4]
    
    MEM_error_KLD_male[j, i]  <-  errors_male[2]
    MEM_error_JS1_male[j, i]  <-  errors_male[3]
    MEM_error_JS2_male[j, i]  <-  errors_male[4]
    
  }
  
  
}

All_errors_MEM<-list( MEM_error_KLD_female,MEM_error_KLD_male,
                         MEM_error_JS1_female,MEM_error_JS1_male,
                         MEM_error_JS2_female,MEM_error_JS2_male)



# saveRDS(All_errors_MEM,paste0(dirl.p,"All_errors_MEM.rds"))




