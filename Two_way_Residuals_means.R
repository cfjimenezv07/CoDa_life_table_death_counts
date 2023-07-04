###############################################################################
#Residuals obtained based on means
###############################################################################
#' Dependencies from other functions
source("FANOVA.R")
#' Parameters for the  Two_way_Residuals_means
#' @param  n.  It represents the total number of functional curves.
#' @param  p.  It represents the grid size.
#' @param  data_pop1. It's a p by n matrix
#' @param  data_pop2. It's a p by n matrix
#' @param  n_prefectures. The number of prefectures, states or departments. 
#' @param  n_year. Number of years considered.
#' @param  n_age. Number of ages considered in each year.
Two_way_Residuals_means<-function(data_pop1,data_pop2
                                  ,year=1959:2020,age= 0:100,n_prefectures=51,n_populations=2){
  FANOVA_means <- FANOVA(data_pop1,data_pop2,year=1959:2020,age= 0:100,n_prefectures=51,n_populations=2)
  # FANOVA_means is the result after running the function FANOVA
  Two_FGE=FANOVA_means$FGE_mean
  Two_FRE=FANOVA_means$FRE_mean
  Two_FCE=FANOVA_means$FCE_mean
  #Number of years considered in each population
  n_year = length(year)
  #Number of ages considered in each year
  n_age = length(age)
  
  
  all_male=t(data_pop1)
  all_female=t(data_pop2)
  
  residuals_b1<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  residuals_b2<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  for (i in 1:nrow(residuals_b1)) {
    residuals_b1[i,]=all_male[i,]-Two_FGE
    residuals_b2[i,]=all_female[i,]-Two_FGE
  }
  
  #remove the column
  residuals_b1c<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  residuals_b2c<-matrix(0,nrow=n_prefectures*n_year,ncol = n_age)
  for (i in 1:nrow(residuals_b1c)) {
    residuals_b1c[i,]=residuals_b1[i,]-Two_FCE[1,]
    residuals_b2c[i,]=residuals_b2[i,]-Two_FCE[2,]
  }
  
  #Remove the row
  residuals_b1r<-matrix(0,nrow=(n_prefectures*n_year),ncol = n_age)
  residuals_b2r<-matrix(0,nrow=(n_prefectures*n_year),ncol = n_age)
  for (j in 1:n_prefectures) {
    residuals_b1r[(n_year*j-(n_year-1)):(n_year*j),]=residuals_b1c[(n_year*j-(n_year-1)):(n_year*j),]-t(replicate(n_year,Two_FRE[j,]))
    residuals_b2r[(n_year*j-(n_year-1)):(n_year*j),]=residuals_b2c[(n_year*j-(n_year-1)):(n_year*j),]-t(replicate(n_year,Two_FRE[j,]))
  }
  
  #Proof Residuals
  data_c1<-matrix(0,nrow=n_prefectures,ncol=n_age)
  data_c2<-matrix(0,nrow=n_prefectures,ncol=n_age)
  for (j in 1:n_prefectures) {
    data_c1[j,]=Two_FGE+Two_FRE[j,]+Two_FCE[1,]
    data_c2[j,]=Two_FGE+Two_FRE[j,]+Two_FCE[2,]
  }
  
  fixed_data_1<-matrix(0,nrow=(n_year*n_prefectures),ncol=n_age)
  fixed_data_2<-matrix(0,nrow=(n_year*n_prefectures),ncol=n_age)
  for (i in 1:n_prefectures) {
    fixed_data_1[(n_year*i-(n_year-1)):(n_year*i),]=do.call(rbind, replicate(n_year, data_c1[i,], simplify=FALSE))
    fixed_data_2[(n_year*i-(n_year-1)):(n_year*i),]=do.call(rbind, replicate(n_year, data_c2[i,], simplify=FALSE))
  }
  FD<-cbind(fixed_data_1,fixed_data_2)
  Recovered_data_1<-fixed_data_1+residuals_b1r
  Recovered_data_2<-fixed_data_2+residuals_b2r
  RD<-cbind(Recovered_data_1,Recovered_data_2)
  
  
  
  #reconstruction proof
  R1=all(round(all_male,4)==round(Recovered_data_1,4))
  R2=all(round(all_female,4)==round(Recovered_data_2,4))
  R<-c(R1,R2)
  
  return(list(residuals1_mean= residuals_b1r,residuals2_mean=residuals_b2r,rd_mean=RD,R_mean=R,Fixed_comp_mean=FD))
  
}

#' Return from Two_way_Residuals_means
#' 
#' @return Return a list containing: residuals from population 1, residuals from population 2, 
#' a logic vector RD, the residuals stacked by columns with both populations, and the fixed components.
#' Everything computed based on means. 
#' @return residuals1_mean, a matrix with  dimension n by p.
#' @return residuals2_mean, a matrix with  dimension n by p.
#' @return rd_mean, a two dimension logic vector that proves that the decomposition sum up to the data.
#' @return R_mean, a matrix of dimension as n by 2p. This represent the time-varying component in the decomposition. 
#' @return Fixed_comp_mean, a matrix of dimension as n by 2p. This represent the deterministic component in the decomposition. 
