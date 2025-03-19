# Uniform Prediction bands

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
# Set a working directory
setwd("/Users/cristianjimenez/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@kaust.edu.sa/My Drive/Spring_2024/STAT_397/PhD_Project_5/Codes_paper2")
##########################
# set a images directory
##########################
dir.l<-"~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@kaust.edu.sa/My Drive/Spring_2024/STAT_397/PhD_Project_5/Plots_paper/"
##########################
# set a results directory
##########################
dirl.p<-"~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@kaust.edu.sa/My Drive/Spring_2024/STAT_397/PhD_Project_5/Resuts_final/USA/"
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

source("Aux_Uniform_prediction_band.R")
Emp_cov_male_80_means_EVR_1lrc<-array(NA,dim = c(n_states,fh=10,2))
Emp_cov_male_95_means_EVR_1lrc<-array(NA,dim = c(n_states,fh=10,2))
Emp_cov_female_80_means_EVR_1lrc<-array(NA,dim = c(n_states,fh=10,2))
Emp_cov_female_95_means_EVR_1lrc<-array(NA,dim = c(n_states,fh=10,2))

Emp_cov_male_80_means_K_1lrc<-array(NA,dim = c(n_states,fh=10,2))
Emp_cov_male_95_means_K_1lrc<-array(NA,dim = c(n_states,fh=10,2))
Emp_cov_female_80_means_K_1lrc<-array(NA,dim = c(n_states,fh=10,2))
Emp_cov_female_95_means_K_1lrc<-array(NA,dim = c(n_states,fh=10,2))


n_states=dataset_entries[[5]]

n_states <- 3
for (i in 1:n_states) {
  
  male_80_means_EVR_1lrc <- coverage_CoDA_uniform(dat=t(male_prefecture_res_means[[i]]),fixed_comp=t(male_prefecture_fixed_means[[i]]),alpha_transf=all_alpha_male[,i],
                                           sample_number=n_year, fh=10, B=1000, level=80, fmethod ="arima",no_core=detectCores()-2,ncomp_selection = "EVR")
  Emp_cov_male_80_means_EVR_1lrc[i,,]<-male_80_means_EVR_1lrc$result
  
  male_80_means_K_1lrc <- coverage_CoDA_uniform(dat=t(male_prefecture_res_means[[i]]),fixed_comp=t(male_prefecture_fixed_means[[i]]),alpha_transf=all_alpha_male[,i],
                                         sample_number=n_year, fh=10, B=1000, level=80, fmethod ="arima",no_core=detectCores()-2,ncomp_selection = "fixed")
  Emp_cov_male_80_means_K_1lrc[i,,]<-male_80_means_K_1lrc$result
  print(i)
}

for (i in 1:n_states) {
  
  male_95_means_EVR_1lrc <- coverage_CoDA_uniform(dat=t(male_prefecture_res_means[[i]]),fixed_comp=t(male_prefecture_fixed_means[[i]]),alpha_transf=all_alpha_male[,i],
                                           sample_number=n_year, fh=10, B=1000, level=95, fmethod ="arima",no_core=detectCores()-2,ncomp_selection = "EVR")
  Emp_cov_male_95_means_EVR_1lrc[i,,]<-male_95_means_EVR_1lrc$result
  
  male_95_means_K_1lrc <- coverage_CoDA_uniform(dat=t(male_prefecture_res_means[[i]]),fixed_comp=t(male_prefecture_fixed_means[[i]]),alpha_transf=all_alpha_male[,i],
                                         sample_number=n_year, fh=10, B=1000, level=95, fmethod ="arima",no_core=detectCores()-2,ncomp_selection = "fixed")
  Emp_cov_male_95_means_K_1lrc[i,,]<-male_95_means_K_1lrc$result 
}

for (i in 1:n_states) {
  
  female_80_means_EVR_1lrc <- coverage_CoDA_uniform(dat=t(female_prefecture_res_means[[i]]),fixed_comp=t(female_prefecture_fixed_means[[i]]),alpha_transf=all_alpha_female[,i],
                                             sample_number=n_year, fh=10, B=1000, level=80, fmethod ="arima",no_core=detectCores()-2,ncomp_selection = "EVR")
  Emp_cov_female_80_means_EVR_1lrc[i,,]<-female_80_means_EVR_1lrc$result
  
  female_80_means_K_1lrc <- coverage_CoDA(dat=t(female_prefecture_res_means[[i]]),fixed_comp=t(female_prefecture_fixed_means[[i]]),alpha_transf=all_alpha_female[,i],
                                          sample_number=n_year, fh=10, B=1000, level=80, fmethod ="arima",no_core=detectCores()-2,ncomp_selection = "fixed")
  Emp_cov_female_80_means_K_1lrc[i,,]<-female_80_means_K_1lrc$result
}

for (i in 1:n_states) {
  
  female_95_means_EVR_1lrc <- coverage_CoDA_uniform(dat=t(female_prefecture_res_means[[i]]),fixed_comp=t(female_prefecture_fixed_means[[i]]),alpha_transf=all_alpha_female[,i],
                                             sample_number=n_year, fh=10, B=1000, level=95, fmethod ="arima",no_core=detectCores()-2,ncomp_selection = "EVR")
  Emp_cov_female_95_means_EVR_1lrc[i,,]<-female_95_means_EVR_1lrc$result
  
  female_95_means_K_1lrc <- coverage_CoDA_uniform(dat=t(female_prefecture_res_means[[i]]),fixed_comp=t(female_prefecture_fixed_means[[i]]),alpha_transf=all_alpha_female[,i],
                                           sample_number=n_year, fh=10, B=1000, level=95, fmethod ="arima",no_core=detectCores()-2,ncomp_selection = "fixed")
  Emp_cov_female_95_means_K_1lrc[i,,]<-female_95_means_K_1lrc$result
  
}

