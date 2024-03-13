# Alternative methods of GSY19 and TNH23 for CoDA
# Replicate Tavakoli et. al. 2023 (THN) and Gao et al. 2019 (GSY) methods for our datasets. 
# These Rcodes are public posted in https://zenodo.org/record/7408999 from TNH23, it includes GSY19.
# We adapted to our paper cases.




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


#####################################################
# 4. Apply the log-ratio transformation
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

all_unconstrained_male     <-t(list.cbind(Unconstrained_male))
all_unconstrained_female   <-t(list.cbind(Unconstrained_female))
all_alpha_male   <- list.cbind(alpha_male)
all_alpha_female <- list.cbind(alpha_female)


USA_male   <- Unconstrained_male 
USA_female <- Unconstrained_female


################################################################################
# auxiliary functions for the computation of the TNH and GSY methods
################################################################################
source('./MITS_class.R') 
source("./aux_TNH.R")


################################################################################
# For the data set: USA female
################################################################################

order <- 3
r <- 3
p <- 9
h_max <- 10 # forecast horizon
years <- sapply(colnames(USA_female[[1]]),toString)
tmp <- init_data_new(data=USA_female,p,Names=names_states)
data_female <- tmp[[1]]  # datasets for Gao2019
mixed_fts_female <- tmp[[2]]
# basis splines
basis <- create.bspline.basis(c(0, 1), nbasis = p, norder = 4)
args <- seq(0, 1, length=n_age)


model_female <- list()
training_set_female <- list()
# Training model for Gao 2019 approach
for(t in 1:h_max){ 
  training_set_female[[years[n_year+t-h_max]]] <- trim_data_new(data_female,t)[[1]]
  model_female[[years[n_year+t-h_max]]] <- hdfpca(training_set_female[[years[n_year+t-h_max]]], order, q = sqrt(dim(training_set_female[[years[n_year+t-h_max]]][[1]])[2]),r) 
}


predictions_GAO_female <- list()
predictions_TNH_female <- list()
for(h in 1:h_max){
  predictions_GAO_female[[h]] <- array(0,dim=c(n_states,h_max-h+1,n_age))
  predictions_TNH_female[[h]] <- array(0,dim=c(n_states,h_max-h+1,n_age))
}

r_hat_female <- list()
cp_female <- list()
for(t in 1:h_max){
  cp_female[[years[n_year+t-h_max]]] <- mixed_fts_female$copy()
  cp_female[[years[n_year+t-h_max]]]$trimSampleSize(n_year+t-h_max)
  r_hat_female[[years[n_year+t-h_max]]] <-   est_r_abc(cp_female[[years[n_year+t-h_max]]],h_max)
  print(paste('estimated r for year ',years[n_year+t-h_max],': ',r_hat_female[[years[n_year+t-h_max]]]))
}


for(h in 1:h_max){
  for (t in 1:(h_max - h + 1)) {
    init_year <- n_year - (h_max-h+1)
    predictions_TNH_female[[h]][,t,] <- predict_NTH_rknown_v2(cp_female[[years[init_year + t]]],h,p,basis,args,r_hat_female[[years[init_year + t]]])
  }
}

for(j in 1:h_max){ 
  tmp <- forecast.hdfpca(model_female[[years[n_year+j-h_max]]],h=h_max - j + 1,B=100)
  pred<- tmp[[1]]
    for(i in 1:n_states){
      pred_state<-pred[[i]] 
      predictions_GAO_female[[j]][i,,] <- t(pred_state)
    }
}

# saveRDS(predictions_GAO_female,paste0(dirl.p,"predictions_GAO_female.rds"))
# saveRDS(predictions_TNH_female,paste0(dirl.p,"predictions_TNH_female.rds"))


forecast_female_TNH <- list()
forecast_female_GSY <- list()
h_max <- 10
for (i in 1:n_states) {
  forecast_female_TNH[[i]] <- array(0, c(n_age, h_max, h_max))
  forecast_female_GSY[[i]] <- array(0, c(n_age, h_max, h_max))
  for (j in 1:h_max) {
    aux <- matrix(0, ncol = j-1, nrow = n_age)
    tnh_aux <- as.matrix(predictions_TNH_female[[j]][i,,]); if (j < h_max) tnh_aux <- t(tnh_aux)
    gao_aux <- as.matrix(predictions_GAO_female[[j]][i,,]); if (j < h_max) gao_aux <- t(gao_aux)
    forecast_female_TNH[[i]][,j,] <- cbind(tnh_aux, aux)
    forecast_female_GSY[[i]][,,j] <- cbind(gao_aux, aux)
  }
}



# saveRDS(forecast_female_TNH,paste0(dirl.p,"forecast_female_TNH.rds"))
# saveRDS(forecast_female_GSY,paste0(dirl.p,"forecast_female_GSY.rds"))

# GSY
h_max <-10
GAO_error_KLD_female<-matrix(0,nrow=h_max,ncol = n_states)
GAO_error_JS1_female<-matrix(0,nrow=h_max,ncol = n_states)
GAO_error_JS2_female<-matrix(0,nrow=h_max,ncol = n_states)
# TNH
TNH_error_KLD_female<-matrix(0,nrow=h_max,ncol = n_states)
TNH_error_JS1_female<-matrix(0,nrow=h_max,ncol = n_states)
TNH_error_JS2_female<-matrix(0,nrow=h_max,ncol = n_states)



Forc_transformed_female_GSY <- list()
Forc_transformed_female_TNH <- list()
for (i in 1:n_states) {

  alpha_x  <- all_alpha_female[,i]
  dat_female <- All_USA_female_qx[[i]] 

  
  forecast_female_TNH_1 <-forecast_female_TNH[[i]]
  forecast_female_GSY_1 <-forecast_female_GSY[[i]]
  # transform forecast
  Forc_transformed_female_TNH_1   <- lapply(1:h_max, transform_back , forecast_female_TNH_1, h_max,alpha_x,n_year,n_age,age)
  Forc_transformed_female_GSY_1   <- lapply(1:h_max, transform_back , forecast_female_GSY_1, h_max,alpha_x,n_year,n_age,age)
  
  Forc_transformed_female_TNH[[i]]   <- build.tensor(l=Forc_transformed_female_TNH_1, nrows=n_age, ncols=h_max,name="Forc_transformed")
  Forc_transformed_female_GSY[[i]]   <- build.tensor(l=Forc_transformed_female_GSY_1, nrows=n_age, ncols=h_max,name="Forc_transformed")
  
  for (j in 1:h_max) {
    
    
    errors_female_TNH  <- rowMeans(sapply( 1:(h_max-j+1),pairwise_errors, Forc_transformed_female_TNH[[i]][,j,],dat_female[,(n_year-h_max+j):(n_year)])) 
    errors_female_GSY  <- rowMeans(sapply( 1:(h_max-j+1),pairwise_errors, Forc_transformed_female_GSY[[i]][,j,],dat_female[,(n_year-h_max+j):(n_year)])) 
 
    
    
    TNH_error_KLD_female[j, i]  <-  errors_female_TNH[1]
    TNH_error_JS1_female[j, i]  <-  errors_female_TNH[2]
    TNH_error_JS2_female[j, i]  <-  errors_female_TNH[3]
    
    GAO_error_KLD_female[j, i]  <-  errors_female_GSY[1]
    GAO_error_JS1_female[j, i]  <-  errors_female_GSY[2]
    GAO_error_JS2_female[j, i]  <-  errors_female_GSY[3]
    
  }

}

All_errors_female<-list( TNH_error_KLD_female,GAO_error_KLD_female,
                         TNH_error_JS1_female,GAO_error_JS1_female,
                         TNH_error_JS2_female,GAO_error_JS2_female)



# saveRDS(All_errors_female,paste0(dirl.p,"All_errors_female_TNH_GSY.rds"))

################################################################################
# For the dataset: USA_male
################################################################################
tmp <- init_data_new(data=USA_male,p,Names=names_states)
data_male <- tmp[[1]]  # datasets for Gao2019
mixed_fts_male <- tmp[[2]]
# basis splines
basis <- create.bspline.basis(c(0, 1), nbasis = p, norder = 4)
args <- seq(0, 1, length=n_age)

model_male <- list()
training_set_male <- list()
# Training model for Gao 2019 approach
for(t in 1:h_max){ 
  training_set_male[[years[n_year+t-h_max]]] <- trim_data_new(data_male,t)[[1]]
  model_male[[years[n_year+t-h_max]]] <- hdfpca(training_set_male[[years[n_year+t-h_max]]], order, q = sqrt(dim(training_set_male[[years[n_year+t-h_max]]][[1]])[2]),r) 
}


predictions_GAO_male <- list()
predictions_TNH_male <- list()
for(h in 1:h_max){
  predictions_GAO_male[[h]] <- array(0,dim=c(n_states,h_max-h+1,n_age))
  predictions_TNH_male[[h]] <- array(0,dim=c(n_states,h_max-h+1,n_age))
}

r_hat_male <- list()
cp_male <- list()
for(t in 1:h_max){
  cp_male[[years[n_year+t-h_max]]] <- mixed_fts_male$copy()
  cp_male[[years[n_year+t-h_max]]]$trimSampleSize(n_year+t-h_max)
  r_hat_male[[years[n_year+t-h_max]]] <-   est_r_abc(cp_male[[years[n_year+t-h_max]]],h_max)
  print(paste('estimated r for year ',years[n_year+t-h_max],': ',r_hat_male[[years[n_year+t-h_max]]]))
}


for(h in 1:h_max){
  for (t in 1:(h_max - h + 1)) {
    init_year <- n_year - (h_max-h+1)
    predictions_TNH_male[[h]][,t,] <- predict_NTH_rknown_v2(cp_male[[years[init_year + t]]],h,p,basis,args,r_hat_male[[years[init_year + t]]])
  }
}

for(j in 1:h_max){ 
  tmp <- forecast.hdfpca(model_male[[years[n_year+j-h_max]]],h=h_max - j + 1,B=100)
  pred<- tmp[[1]]
  for(i in 1:n_states){
    pred_state<-pred[[i]] 
    predictions_GAO_male[[j]][i,,] <- t(pred_state)
  }
}

# saveRDS(predictions_GAO_male,paste0(dirl.p,"predictions_GAO_male.rds"))
# saveRDS(predictions_TNH_male,paste0(dirl.p,"predictions_TNH_male.rds"))



forecast_male_TNH <- list()
forecast_male_GSY <- list()

for (i in 1:n_states) {
  forecast_male_TNH[[i]] <- array(0, c(n_age, h_max, h_max))
  forecast_male_GSY[[i]] <- array(0, c(n_age, h_max, h_max))
  for (j in 1:h_max) {
    aux <- matrix(0, ncol = j-1, nrow = n_age)
    tnh_aux <- as.matrix(predictions_TNH_male[[j]][i,,]); if (j < h_max) tnh_aux <- t(tnh_aux)
    gao_aux <- as.matrix(predictions_GAO_male[[j]][i,,]); if (j < h_max) gao_aux <- t(gao_aux)
    forecast_male_TNH[[i]][,j,] <- cbind(tnh_aux, aux)
    forecast_male_GSY[[i]][,,j] <- cbind(gao_aux, aux)
  }
}

# saveRDS(forecast_male_TNH,paste0(dirl.p,"forecast_male_TNH.rds"))
# saveRDS(forecast_male_GSY,paste0(dirl.p,"forecast_male_GSY.rds"))

# GSY

GAO_error_KLD_male<-matrix(0,nrow=h_max,ncol = n_states)
GAO_error_JS1_male<-matrix(0,nrow=h_max,ncol = n_states)
GAO_error_JS2_male<-matrix(0,nrow=h_max,ncol = n_states)
# TNH
TNH_error_KLD_male<-matrix(0,nrow=h_max,ncol = n_states)
TNH_error_JS1_male<-matrix(0,nrow=h_max,ncol = n_states)
TNH_error_JS2_male<-matrix(0,nrow=h_max,ncol = n_states)



Forc_transformed_male_GSY <- list()
Forc_transformed_male_TNH <- list()
for (i in 1:n_states) {
  
  alpha_x  <- all_alpha_male[,i]
  dat_male <- All_USA_male_qx[[i]] 
  
  forecast_male_TNH_1 <-forecast_male_TNH[[i]]
  forecast_male_GSY_1 <-forecast_male_GSY[[i]]
  # transform forecast
  Forc_transformed_male_TNH_1   <- lapply(1:h_max, transform_back , forecast_male_TNH_1, h_max,alpha_x,n_year,n_age,age)
  Forc_transformed_male_GSY_1   <- lapply(1:h_max, transform_back , forecast_male_GSY_1, h_max,alpha_x,n_year,n_age,age)
  
  Forc_transformed_male_TNH[[i]]   <- build.tensor(l=Forc_transformed_male_TNH_1, nrows=n_age, ncols=h_max,name="Forc_transformed")
  Forc_transformed_male_GSY[[i]]   <- build.tensor(l=Forc_transformed_male_GSY_1, nrows=n_age, ncols=h_max,name="Forc_transformed")
  
  for (j in 1:h_max) {
    
    
    errors_male_TNH  <- rowMeans(sapply( 1:(h_max-j+1),pairwise_errors, Forc_transformed_male_TNH[[i]][,j,],dat_male[,(n_year-h_max+j):(n_year)])) 
    errors_male_GSY  <- rowMeans(sapply( 1:(h_max-j+1),pairwise_errors, Forc_transformed_male_GSY[[i]][,j,],dat_male[,(n_year-h_max+j):(n_year)])) 
    
    
    
    TNH_error_KLD_male[j, i]  <-  errors_male_TNH[1]
    TNH_error_JS1_male[j, i]  <-  errors_male_TNH[2]
    TNH_error_JS2_male[j, i]  <-  errors_male_TNH[3]
    
    GAO_error_KLD_male[j, i]  <-  errors_male_GSY[1]
    GAO_error_JS1_male[j, i]  <-  errors_male_GSY[2]
    GAO_error_JS2_male[j, i]  <-  errors_male_GSY[3]
    
  }
  
  
}

All_errors_male<-list( TNH_error_KLD_male,GAO_error_KLD_male,
                         TNH_error_JS1_male,GAO_error_JS1_male,
                         TNH_error_JS2_male,GAO_error_JS2_male)



saveRDS(All_errors_male,paste0(dirl.p,"All_errors_male_TNH_GSY.rds"))

v1 <- round(cbind(apply(All_errors_female[[1]],1,mean),apply(All_errors_female[[5]],1,mean),apply(All_errors_male[[1]],1,mean),
            apply(All_errors_male[[5]],1,mean),
            +             apply(All_errors_female[[2]],1,mean),apply(All_errors_female[[6]],1,mean),
            apply(All_errors_male[[2]],1,mean),apply(All_errors_male[[6]],1,mean))*100,2)
