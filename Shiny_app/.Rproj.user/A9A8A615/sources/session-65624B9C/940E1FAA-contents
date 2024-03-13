#Pvdata
library('fda')
library('fda.usc')
library('signal')
library('pracma')
library(roahd)
library(readxl)
library(fdaoutlier)
library(RandomFields)
library(ggplot2)
library(reshape2)
library(grid)
library(scatterplot3d)
library(pvclust)
library(MASS)
library(fda.usc)
library(ddalpha)
library(fdaoutlier)
`%!in%`=Negate(`%in%`)
#Outliers from fbplot
outliers_fbplot<-function(data){
  msp= fda::fbplot((data),plot=F)
  outliers= msp$outpoint
  return(outliers)
}
# Outliers from MS plot
outliers_msplot<-function(data){
  msp= fdaoutlier::functional_boxplot(t(data),
                                      depth_method='dirout')
  outliers= msp$outliers
  return(outliers)
}
# function to obtain the median from Fbplot
Get_median<-function(mod){
  dd=fda::fbplot(mod,plot=F)$depth
  med=which(dd==max(dd))
  return(med)
}
combinat=function(n,p){if (n<p){combinat=0}else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}}
PV.hourly.Data=read.csv("~/Google Drive/Spring 2022/STAT 397/PhD project 1/Rcodes/PV hourly Data - 2015-2020.csv",header=T,sep = ',')
dir<-"/Volumes/GoogleDrive/My Drive/Spring 2022/STAT 397/PhD project 1/Rcodes"
PV.hourly.Data=PV.hourly.Data[,-c(5,6)]
# PV_SI=cor(PV.hourly.Data[,2],PV.hourly.Data[,3])
# PV_Sun=cor(PV.hourly.Data[,2],PV.hourly.Data[,4])
# SI_Sun=cor(PV.hourly.Data[,3],PV.hourly.Data[,4])

M=10 #Daily sun hours
N=2192 
PV.data=matrix(0,nrow=M,ncol =N)
G.data=matrix(0,nrow=M,ncol = N)
Sun.data=matrix(0,nrow=M,ncol = N)
for (i in 1:N) {
  PV.data[,i]=PV.hourly.Data[c((24*i-18):(24*i-9)),2]
  G.data[,i]=PV.hourly.Data[c((24*i-18):(24*i-9)),3]
  Sun.data[,i]=PV.hourly.Data[c((24*i-18):(24*i-9)),4]
}



PV.data_1=PV.data[,1:365]
PV.data_2=PV.data[,366:731]
PV.data_3=PV.data[,732:1096]
PV.data_4=PV.data[,1097:1461]
PV.data_5=PV.data[,1462:1826]
PV.data_6=PV.data[,1827:2192]


#From PV power

fm_2019<-cbind(PV.data_5[,15:166])
fm_2020<-cbind(PV.data_6[,15:166])
winter_PV<-list(fm_2019,fm_2020)


#Outliers detection procedure
par(mfrow=c(1,1))
outliers_PV<-list()
outliers_magPV<-list()
for (h in 1:2) {
  data=winter_PV[[h]]
  outl_fbplot=outliers_fbplot(data)
  outliers_magPV[[h]]<-outl_fbplot
  med=Get_median(data)
  M=dim(data)[1]
  N=dim(data)[2]
  t=seq(9,16,len=M)
  plot(NaN,xlim=c(8,17),ylim=c(min(data),max(data)),ylab='PV system [W]',xlab='Time [h] ')
  for (i in 1:N) {
    if(i%!in%outl_fbplot){
      lines(t,data[,i],col='azure3')
    }
  }
  
  for (i in 1:N) {
    if(i%in%outl_fbplot){
      lines(t,data[,i],col='darkturquoise')
    }
  }
  for (i in 1:N) {
    if(i%in%med){
      lines(t,data[,i],col=9,lwd=2)
    }
  }
  data=data[,-outl_fbplot]
  M=dim(data)[1]
  N=dim(data)[2]
  t=seq(0,1,len=M)
  med=Get_median(data)
  Pw.rank=as.matrix(t(apply(data,1,rank))) 
  abo=N-Pw.rank
  bel=Pw.rank-1
  Pw.depth=as.matrix((abo*bel)/combinat(N,2)+(N-1)/combinat(N,2)) 
  Pair.Pw.depth=matrix(NaN,nrow =(M-1),ncol = 2*N)
  for (i in 1:dim(data)[2]) {
    Pair.Pw.depth[,(2*i-1)]= Pw.depth[1:(M-1),i]
    Pair.Pw.depth[,(2*i)]= Pw.depth[2:M,i]
  }
  prod.var=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2) )
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    prod.var[i,]=var(Pair.Pw.depth[,(2*i-1)])*var(Pair.Pw.depth[,(2*i)])
  }
  #Bivariate Distribution of Non-Outlying curves
  corr=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2))
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    if( prod.var[i,]!=0){
      corr[i,]=cor(Pair.Pw.depth[,(2*i-1)],Pair.Pw.depth[,(2*i)])
    }else{
      corr[i,]=0
    }
  }
  #Boxplot of sample correlation
  a=boxplot.stats(corr)[[1]][1]
  outl=as.vector(which(corr<=a))
  outliers_PV[[h]]<-outl
  #data
  t=seq(9,16,len=M)
  plot(NaN,xlim=c(8,17),ylim=c(min(data),max(data)),ylab='PV system [W]',xlab='Time [h] ')
  for (i in 1:N) {
    if(i%!in%outl){
      lines(t,data[,i],col='azure3')
    }
  }
  for (i in 1:N) {
    if(i%in%outl){
      lines(t,data[,i],col='darkmagenta',lwd=2)
    }
  }
  for (i in 1:N) {
    if(i%in%med){
      lines(t,data[,i],col=9,lwd=2)
    }
  }
  #Bivariate distribution

}

################################################################################
#From Global solar irrandiance,
SI.data_5=G.data[,1462:1826]
SI.data_6=G.data[,1827:2192]

fm_SI_2019<-cbind(SI.data_5[,15:166])
fm_SI_2020<-cbind(SI.data_6[,15:166])
winter_SI<-list(fm_SI_2019,fm_SI_2020)


#from fbplot
outliers_SI<-list()
#par(mfrow=c(1,1))
outliers_mag_SI<-list()
for (h in 1:2) {
  data=winter_SI[[h]]
  outl_fbplot=outliers_fbplot(data)
  outliers_mag_SI[[h]]<- outl_fbplot
  med=Get_median(data)
  M=dim(data)[1]
  N=dim(data)[2]
  t=seq(9,16,len=M)
  plot(NaN,xlim=c(8,17),ylim=c(min(data),max(data)),ylab='Global solar irradiance [W/m^2]',xlab='Time [h] ')
  for (i in 1:N) {
    if(i%!in%outl_fbplot){
      lines(t,data[,i],col='azure3')
    }
  }
  
  for (i in 1:N) {
    if(i%in%outl_fbplot){
      lines(t,data[,i],col='darkturquoise')
    }
  }
  for (i in 1:N) {
    if(i%in%med){
      lines(t,data[,i],col=9,lwd=2)
    }
  }
  data=data[,-outl_fbplot]
  M=dim(data)[1]
  N=dim(data)[2]
  t=seq(0,1,len=M)
  med=Get_median(data)
  Pw.rank=as.matrix(t(apply(data,1,rank))) 
  abo=N-Pw.rank
  bel=Pw.rank-1
  Pw.depth=as.matrix((abo*bel)/combinat(N,2)+(N-1)/combinat(N,2)) 
  Pair.Pw.depth=matrix(NaN,nrow =(M-1),ncol = 2*N)
  for (i in 1:dim(data)[2]) {
    Pair.Pw.depth[,(2*i-1)]= Pw.depth[1:(M-1),i]
    Pair.Pw.depth[,(2*i)]= Pw.depth[2:M,i]
  }
  prod.var=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2) )
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    prod.var[i,]=var(Pair.Pw.depth[,(2*i-1)])*var(Pair.Pw.depth[,(2*i)])
  }
  #Bivariate Distribution of Non-Outlying curves
  corr=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2))
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    if( prod.var[i,]!=0){
      corr[i,]=cor(Pair.Pw.depth[,(2*i-1)],Pair.Pw.depth[,(2*i)])
    }else{
      corr[i,]=0
    }
  }
  #Boxplot of sample correlation
  a=boxplot.stats(corr)[[1]][1]
  outl=as.vector(which(corr<=a))
  outliers_SI[[h]]<-outl
  #data
  t=seq(9,16,len=M)
  plot(NaN,xlim=c(8,17),ylim=c(min(data),max(data)),ylab='Global solar irradiance [W/m^2]',xlab='Time [h] ')
  for (i in 1:N) {
    if(i%!in%outl){
      lines(t,data[,i],col='azure3')
    }
  }
  for (i in 1:N) {
    if(i%in%outl){
      lines(t,data[,i],col='darkmagenta',lwd=2)
    }
  }
  for (i in 1:N) {
    if(i%in%med){
      lines(t,data[,i],col=9,lwd=2)
    }
  }
}


nn<-dim(fm_SI_2020)[2]
pdo1_PV<-(length(outliers_PV[[1]])/(nn-length(outliers_magPV[[1]])))*100
pdo2_PV<-(length(outliers_PV[[2]])/(nn-length(outliers_magPV[[2]])))*100
pdo1_SI<-(length(outliers_SI[[1]])/(nn-length(outliers_mag_SI[[1]])))*100
pdo2_SI<-(length(outliers_SI[[2]])/(nn-length(outliers_mag_SI[[2]])))*100


################################################################################

#Bivariate Functional data PV-Solar Irradiance

#PV data
fm_2019<-cbind(PV.data_5[,15:166])
fm_2020<-cbind(PV.data_6[,15:166])


#Solar Irradiance
fm_SI_2019<-cbind(SI.data_5[,15:166])
fm_SI_2020<-cbind(SI.data_6[,15:166])


#magnitude outliers
Magnitude_outliers_both_1<-unique(outliers_mag_SI[[1]],outliers_magPV[[1]])
Magnitude_outliers_both_2<-unique(outliers_mag_SI[[2]],outliers_magPV[[2]])
PV<-list(fm_2019[,-Magnitude_outliers_both_1],fm_2020[,-Magnitude_outliers_both_2])
G.data<-list(fm_SI_2019[,-Magnitude_outliers_both_1],fm_SI_2020[,-Magnitude_outliers_both_2])


#par(mfrow=c(1,1))
Biv_data=list()
med2=list()
p=dim(fm_SI_2019)[1]
n=dim(fm_SI_2019)[2]
d=2
for (m in 1:length(PV)) {
  data_PV=PV[[m]]
  data_SI=G.data[[m]]
  #PV-Solar Irradiance
  n=dim(data_PV)[2]
  biv_data = array(0, dim = c(n, p, d))
  biv_data[,,1]=t(as.matrix(data_PV))
  biv_data[,,2]=t(as.matrix(data_SI))
  Biv_data[[m]]=biv_data
  #plot_data(biv_data2)
  med2[[m]]=fdaoutlier::msplot(biv_data,plot = F)$median_curve
} 

#remove magnitude outliers
outliers<-list()
for (h in 1:2) {
  mod=Biv_data[[h]]
  #Computing outliers through the simplicial
  n=dim(mod)[1]
  PSD=matrix(0,nrow=p,ncol = n)
  for (i in 1:n) {
    SD=mdepth.SD(mod[i,,],mod[i,,])
    PSD[,i]=SD$dep
  }
  PSDm<-apply(PSD,2,mean)
  # out_1=mag_outliers_fromSD(mod[,,1],PSDm)
  # out_2=mag_outliers_fromSD(mod[,,2],PSDm)
  
  Pair.Pw.depth=matrix(NaN,nrow =(p-1),ncol = 2*n)
  for (i in 1:dim(PSD)[2]) {
    Pair.Pw.depth[,(2*i-1)]= PSD[1:(p-1),i]
    Pair.Pw.depth[,(2*i)]= PSD[2:p,i]
  }
  
  prod.var=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2) )
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    prod.var[i,]=var(Pair.Pw.depth[,(2*i-1)])*var(Pair.Pw.depth[,(2*i)])
    
  }
  
  #Bivariate Distribution of Non-Outlying curves
  corr=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2))
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    if( prod.var[i,]!=0){
      corr[i,]=cor(Pair.Pw.depth[,(2*i-1)],Pair.Pw.depth[,(2*i)])
    }else{
      corr[i,]=0
    }
  }
  
  a=boxplot.stats(corr)[[1]][1]
  outliers[[h]]=as.vector(which(corr<=a))
  
}

n<-dim(Biv_data[[1]])[1]

pdo1_PVSI<-(length(outliers[[1]])/(n-length(Magnitude_outliers_both_1)))*100
pdo2_PVSI<-(length(outliers[[2]])/(n-length(Magnitude_outliers_both_2)))*100


plot_data = function (data = data){
  x = seq(8, 17) 
  y = data[1, , 2]
  z = data[1, , 1]
  s3d = scatterplot3d(x, y, z,
                      type = "l",
                      pch = 19,
                      color = 'gray86',
                      lab = c(1, 1, 2),
                      angle = 45,
                      mar = c(3,3,3,3.5),
                      xlim = c(7,17),
                      ylim = c(0, 2000),
                      zlim = c(0,250000),
                      xlab = "Time [h]",
                      ylab = "Global solar irradiance [W/m^2]",
                      zlab = "PV power [W]",
                      #main = " Bivariate  Curves",
                      #font.axis = 0.2,
                      cex.lab = 0.7,
                      box = FALSE,
                      axis = T,
                      tick.marks=T,
                      label.tick.marks=T
  )
  for (i in 1:n){
    if(i%!in% o){
      x = seq(8, 17) 
      y = data[i, , 2]
      z = data[i, , 1]
      s3d$points3d(x, 
                   y, 
                   z, 
                   pch = 19, 
                   col = 'azure3',
                   type = "l")
    }
  }
  for (i in 1:n){
    if(i%in% o){
      x = seq(8, 17) 
      y = data[i, , 2]
      z = data[i, , 1]
      s3d$points3d(x, 
                   y, 
                   z, 
                   pch = 19, 
                   col = 'darkmagenta',
                   type = "l",lwd=2)
    }
  }
}

par(mfrow=c(1,1))
m=1
data=Biv_data[[m]]
o=outliers[[m]]
n=dim(data)[1]
plot_data(data)

m=2
data=Biv_data[[m]]
o=outliers[[m]]
n=dim(data)[1]
plot_data(data)

#Plots for the paper
par(mfrow=c(2,2))
################################################################################
#plot_1 Univ PV magnitude
outliers_PV<-list()
outliers_magPV<-list()
par(mfrow=c(2,2))
for (h in 1:1) {
  data=winter_PV[[h]]
  outl_fbplot=outliers_fbplot(data)
  outliers_magPV[[h]]<-outl_fbplot
  med=Get_median(data)
  M=dim(data)[1]
  N=dim(data)[2]
  t=seq(9,16,len=M)
  plot(NaN,xlim=c(8,17),ylim=c(min(data),max(data)),ylab='PV power [W]',xlab='Time [h]')
  for (i in 1:N) {
    if(i%!in%outl_fbplot){
      lines(t,data[,i],col='azure3')
    }
  }
  
  for (i in 1:N) {
    if(i%in%outl_fbplot){
      lines(t,data[,i],col='darkturquoise')
    }
  }
  for (i in 1:N) {
    if(i%in%med){
      lines(t,data[,i],col=9,lwd=2)
    }
  }
  data=data[,-outl_fbplot]
  M=dim(data)[1]
  N=dim(data)[2]
  t=seq(0,1,len=M)
  med=Get_median(data)
  Pw.rank=as.matrix(t(apply(data,1,rank))) 
  abo=N-Pw.rank
  bel=Pw.rank-1
  Pw.depth=as.matrix((abo*bel)/combinat(N,2)+(N-1)/combinat(N,2)) 
  Pair.Pw.depth=matrix(NaN,nrow =(M-1),ncol = 2*N)
  for (i in 1:dim(data)[2]) {
    Pair.Pw.depth[,(2*i-1)]= Pw.depth[1:(M-1),i]
    Pair.Pw.depth[,(2*i)]= Pw.depth[2:M,i]
  }
  prod.var=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2) )
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    prod.var[i,]=var(Pair.Pw.depth[,(2*i-1)])*var(Pair.Pw.depth[,(2*i)])
  }
  #Bivariate Distribution of Non-Outlying curves
  corr=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2))
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    if( prod.var[i,]!=0){
      corr[i,]=cor(Pair.Pw.depth[,(2*i-1)],Pair.Pw.depth[,(2*i)])
    }else{
      corr[i,]=0
    }
  }
  #Boxplot of sample correlation
  a=boxplot.stats(corr)[[1]][1]
  outl=as.vector(which(corr<=a))
  outliers_PV[[h]]<-outl
  #data
  t=seq(9,16,len=M)
  plot(NaN,xlim=c(8,17),ylim=c(min(data),max(data)),ylab='PV power [W]',xlab='Time [h]')
  for (i in 1:N) {
    if(i%!in%outl){
      lines(t,data[,i],col='azure3')
    }
  }
  for (i in 1:N) {
    if(i%in%outl){
      lines(t,data[,i],col='darkmagenta',lwd=2)
    }
  }
  for (i in 1:N) {
    if(i%in%med){
      lines(t,data[,i],col=9,lwd=2)
    }
  }
  #Bivariate distribution
  
  #boxplot of the correlation values
  #boxplot(corr,outcol='6',ylab='SC',pch=20)
}

#2
outliers_SI<-list()
outliers_mag_SI<-list()
for (h in 1:1) {
  data=winter_SI[[h]]
  outl_fbplot=outliers_fbplot(data)
  outliers_mag_SI[[h]]<- outl_fbplot
  med=Get_median(data)
  M=dim(data)[1]
  N=dim(data)[2]
  t=seq(9,16,len=M)
  plot(NaN,xlim=c(8,17),ylim=c(min(data),max(data)),ylab='Global solar irradiance [W/m^2]',xlab='Time [h]')
  for (i in 1:N) {
    if(i%!in%outl_fbplot){
      lines(t,data[,i],col='azure3')
    }
  }
  
  for (i in 1:N) {
    if(i%in%outl_fbplot){
      lines(t,data[,i],col='darkturquoise')
    }
  }
  for (i in 1:N) {
    if(i%in%med){
      lines(t,data[,i],col=9,lwd=2)
    }
  }
  data=data[,-outl_fbplot]
  M=dim(data)[1]
  N=dim(data)[2]
  t=seq(0,1,len=M)
  med=Get_median(data)
  Pw.rank=as.matrix(t(apply(data,1,rank))) 
  abo=N-Pw.rank
  bel=Pw.rank-1
  Pw.depth=as.matrix((abo*bel)/combinat(N,2)+(N-1)/combinat(N,2)) 
  Pair.Pw.depth=matrix(NaN,nrow =(M-1),ncol = 2*N)
  for (i in 1:dim(data)[2]) {
    Pair.Pw.depth[,(2*i-1)]= Pw.depth[1:(M-1),i]
    Pair.Pw.depth[,(2*i)]= Pw.depth[2:M,i]
  }
  prod.var=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2) )
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    prod.var[i,]=var(Pair.Pw.depth[,(2*i-1)])*var(Pair.Pw.depth[,(2*i)])
  }
  #Bivariate Distribution of Non-Outlying curves
  corr=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2))
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    if( prod.var[i,]!=0){
      corr[i,]=cor(Pair.Pw.depth[,(2*i-1)],Pair.Pw.depth[,(2*i)])
    }else{
      corr[i,]=0
    }
  }
  #Boxplot of sample correlation
  a=boxplot.stats(corr)[[1]][1]
  outl=as.vector(which(corr<=a))
  outliers_SI[[h]]<-outl
  #data
  t=seq(9,16,len=M)
  plot(NaN,xlim=c(8,17),ylim=c(min(data),max(data)),ylab='Global solar irradiance [W/m^2]',xlab='Time [h]')
  for (i in 1:N) {
    if(i%!in%outl){
      lines(t,data[,i],col='azure3')
    }
  }
  for (i in 1:N) {
    if(i%in%outl){
      lines(t,data[,i],col='darkmagenta',lwd=2)
    }
  }
  for (i in 1:N) {
    if(i%in%med){
      lines(t,data[,i],col=9,lwd=2)
    }
  }
  #Bivariate distribution
  
  #boxplot of the correlation values
  #boxplot(corr,outcol='6',ylab='SC',pch=20)
}




################################################################################
par(mfrow=c(2,2))
#from fbplot
outliers_PV<-list()
outliers_magPV<-list()
for (h in 2:2) {
  data=winter_PV[[h]]
  outl_fbplot=outliers_fbplot(data)
  outliers_magPV[[h]]<-outl_fbplot
  med=Get_median(data)
  M=dim(data)[1]
  N=dim(data)[2]
  t=seq(9,16,len=M)
  plot(NaN,xlim=c(8,17),ylim=c(min(data),max(data)),ylab='PV power [W]',xlab='Time [h]')
  for (i in 1:N) {
    if(i%!in%outl_fbplot){
      lines(t,data[,i],col='azure3')
    }
  }
  
  for (i in 1:N) {
    if(i%in%outl_fbplot){
      lines(t,data[,i],col='darkturquoise')
    }
  }
  for (i in 1:N) {
    if(i%in%med){
      lines(t,data[,i],col=9,lwd=2)
    }
  }
  data=data[,-outl_fbplot]
  M=dim(data)[1]
  N=dim(data)[2]
  t=seq(0,1,len=M)
  med=Get_median(data)
  Pw.rank=as.matrix(t(apply(data,1,rank))) 
  abo=N-Pw.rank
  bel=Pw.rank-1
  Pw.depth=as.matrix((abo*bel)/combinat(N,2)+(N-1)/combinat(N,2)) 
  Pair.Pw.depth=matrix(NaN,nrow =(M-1),ncol = 2*N)
  for (i in 1:dim(data)[2]) {
    Pair.Pw.depth[,(2*i-1)]= Pw.depth[1:(M-1),i]
    Pair.Pw.depth[,(2*i)]= Pw.depth[2:M,i]
  }
  prod.var=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2) )
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    prod.var[i,]=var(Pair.Pw.depth[,(2*i-1)])*var(Pair.Pw.depth[,(2*i)])
  }
  #Bivariate Distribution of Non-Outlying curves
  corr=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2))
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    if( prod.var[i,]!=0){
      corr[i,]=cor(Pair.Pw.depth[,(2*i-1)],Pair.Pw.depth[,(2*i)])
    }else{
      corr[i,]=0
    }
  }
  #Boxplot of sample correlation
  a=boxplot.stats(corr)[[1]][1]
  outl=as.vector(which(corr<=a))
  outliers_PV[[h]]<-outl
  #data
  t=seq(9,16,len=M)
  plot(NaN,xlim=c(8,17),ylim=c(min(data),max(data)),ylab='PV power [W]',xlab='Time [h]')
  for (i in 1:N) {
    if(i%!in%outl){
      lines(t,data[,i],col='azure3')
    }
  }
  for (i in 1:N) {
    if(i%in%outl){
      lines(t,data[,i],col='darkmagenta',lwd=2)
    }
  }
  for (i in 1:N) {
    if(i%in%med){
      lines(t,data[,i],col=9,lwd=2)
    }
  }
  #Bivariate distribution
  
  #boxplot of the correlation values
  #boxplot(corr,outcol='6',ylab='SC',pch=20)
}

outliers_SI<-list()
outliers_mag_SI<-list()
for (h in 2:2) {
  data=winter_SI[[h]]
  outl_fbplot=outliers_fbplot(data)
  outliers_mag_SI[[h]]<- outl_fbplot
  med=Get_median(data)
  M=dim(data)[1]
  N=dim(data)[2]
  t=seq(9,16,len=M)
  plot(NaN,xlim=c(8,17),ylim=c(min(data),max(data))
       ,ylab='Global solar irradiance [W/m^2]',xlab='Time [h]')
  for (i in 1:N) {
    if(i%!in%outl_fbplot){
      lines(t,data[,i],col='azure3')
    }
  }
  
  for (i in 1:N) {
    if(i%in%outl_fbplot){
      lines(t,data[,i],col='darkturquoise')
    }
  }
  for (i in 1:N) {
    if(i%in%med){
      lines(t,data[,i],col=9,lwd=2)
    }
  }
  data=data[,-outl_fbplot]
  M=dim(data)[1]
  N=dim(data)[2]
  t=seq(0,1,len=M)
  med=Get_median(data)
  Pw.rank=as.matrix(t(apply(data,1,rank))) 
  abo=N-Pw.rank
  bel=Pw.rank-1
  Pw.depth=as.matrix((abo*bel)/combinat(N,2)+(N-1)/combinat(N,2)) 
  Pair.Pw.depth=matrix(NaN,nrow =(M-1),ncol = 2*N)
  for (i in 1:dim(data)[2]) {
    Pair.Pw.depth[,(2*i-1)]= Pw.depth[1:(M-1),i]
    Pair.Pw.depth[,(2*i)]= Pw.depth[2:M,i]
  }
  prod.var=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2) )
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    prod.var[i,]=var(Pair.Pw.depth[,(2*i-1)])*var(Pair.Pw.depth[,(2*i)])
  }
  #Bivariate Distribution of Non-Outlying curves
  corr=matrix(0,ncol = 1,nrow = (dim(Pair.Pw.depth)[2]/2))
  for (i in 1:(dim(Pair.Pw.depth)[2]/2)) {
    if( prod.var[i,]!=0){
      corr[i,]=cor(Pair.Pw.depth[,(2*i-1)],Pair.Pw.depth[,(2*i)])
    }else{
      corr[i,]=0
    }
  }
  #Boxplot of sample correlation
  a=boxplot.stats(corr)[[1]][1]
  outl=as.vector(which(corr<=a))
  outliers_SI[[h]]<-outl
  #data
  t=seq(9,16,len=M)
  plot(NaN,xlim=c(8,17),ylim=c(min(data),max(data)),ylab='Global solar irradiance [W/m^2]',xlab='Time [h]')
  for (i in 1:N) {
    if(i%!in%outl){
      lines(t,data[,i],col='azure3')
    }
  }
  for (i in 1:N) {
    if(i%in%outl){
      lines(t,data[,i],col='darkmagenta',lwd=2)
    }
  }
  for (i in 1:N) {
    if(i%in%med){
      lines(t,data[,i],col=9,lwd=2)
    }
  }
  #Bivariate distribution
  
  #boxplot of the correlation values
  #boxplot(corr,outcol='6',ylab='SC',pch=20)
}



#bivariate
par(mfrow=c(1,2))
for (m in 1:2) {
  data=Biv_data[[m]]
  o=outliers[[m]]
  n=dim(data)[1]
  plot_data(data)
}
