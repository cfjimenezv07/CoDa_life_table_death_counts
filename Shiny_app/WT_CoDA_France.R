# FMP and FM to the constrained densities
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
setwd("~/My Drive/Fall 2023/STAT 397/PhD project 5/Rcodes/")

##########################
# set a images directory
##########################
dir.l<-"~/My Drive/Fall 2023/STAT 397/PhD project 5/Plots_paper/"

################################################################################
# Dataset entries
################################################################################
year = 1968:2021
n_year = length(year)
age = 0:105
n_age = length(age)
n_prefectures=95

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

################################################################################
# Dataset 
################################################################################
name_departments <- readRDS("./name_departments.rds")
All_France_male_qx   <- readRDS("./All_France_male_qx.rds")
All_France_female_qx <- readRDS("./All_France_female_qx.rds")


source("Two_way_FM_Res.R")
France_male_qx <-   t(list.cbind(All_France_male_qx))
France_female_qx <- t(list.cbind(All_France_female_qx))
Y <- cbind(France_male_qx,France_female_qx)

# Compute the functional median polish decomposition.
FMP = ftsa::Two_way_median_polish(Y,n_prefectures = n_prefectures, year = year, age =age, n_populations = n_populations )
##1. The funcional grand effect
FGE = FMP$grand_effect
##2. The funcional row effect
FRE = FMP$row_effect
##3. The funcional column effect
FCE = FMP$col_effect

source("Two_way_FM_Res.R")

# Compute the residuals
Residuals <- Two_way_Residuals(cbind(France_male_qx, France_female_qx),FMP,n_prefectures,n_year,n_age)
Res1=Residuals$residuals1
Res2=Residuals$residuals2
Residuals_<-cbind(Res1,Res2)
# Reconstructed data
RR<-Residuals$rd #Matrix with the original data reconstructed from the FMP decomposition
#  It's the proof of the reconstruction of the residuals. 
Residuals$R #The result should be a vector with two entries TRUE, TRUE.
Fixed_part<-Residuals$Fixed_comp #Fixed components to be added up after forecasting

################################################################################
# Forecasting
################################################################################
source("Point_forecasting.R")

################################################################################
# With a different horizon rolling
################################################################################
forecasted_curves_triangular_France<-list()
for (i in 1:n_prefectures) {
  max_h=10
  pref=(n_year*i-(n_year-1)):(n_year*i)
  forecasted_prefecture<-matrix(0,nrow = n_age*2,ncol=max_h)
  for (k in 1:max_h) {
    pref_k<-pref[k:(n_year-max_h-1+k)]
    t = Sys.time()
    frc=Pref_forecasted_curves(fixed_com=Fixed_part[pref_k,],
                               Residuals_f=Residuals_[pref_k,],
                               est_method = "lrc",
                               fh = 1, PI = NULL, B = 1000)
    forecasted_prefecture[,k]=frc$med_polish_curve_forecast
    print(Sys.time() - t)
    print(k)
  }
  forecasted_curves_triangular_France[[i]]=forecasted_prefecture
  print(i)
}



