# Plots for the Gini coefficients. Generates Figures 2.a and 2.b in the manuscript

# Set the working directory
##########################
# set a working directory
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

savepdf <- function(file, width=16, height=10)
{
  fname <- paste(dir.l,file,".pdf",sep="")
  pdf(fname, width=width/2.54, height=height/2.54,
      pointsize=10)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}

# Gini coefficients
Gini_female_USA <- readRDS("./Gini_female_USA.rds")
Gini_male_USA <- readRDS("./Gini_male_USA.rds")
year_USA = 1959:2020
# year_France = 1968:2021
# USA female
Summary_USA_female     <- matrix(NA,nrow=62,ncol=4)
Summary_USA_female[,1] <- apply(Gini_female_USA, 1,min)
Summary_USA_female[,2] <- apply(Gini_female_USA, 1,median)
Summary_USA_female[,3] <- apply(Gini_female_USA, 1,mean)
Summary_USA_female[,4] <- apply(Gini_female_USA, 1,max)

savepdf("Fig_2a")
matplot(year_USA,Summary_USA_female,type = 'l',lty = c(1,1,1,1),col=c(2,3,4,5),
        main="Female",ylab="Gini coefficient",xlab="Year"
        ,ylim = c(min(Summary_USA_female,Summary_USA_male),max(Summary_USA_female,Summary_USA_male)+0.09))
legend(x = "topright", lty = c(1,1,1,1), cex = 1, 
       col= c(5,4,3,2),text.col = 1, 
       legend=c("max", "mean","median","min"))
dev.off()

# USA male
Summary_USA_male     <- matrix(NA,nrow=62,ncol=4)
Summary_USA_male[,1] <- apply(Gini_male_USA, 1,min)
Summary_USA_male[,2] <- apply(Gini_male_USA, 1,median)
Summary_USA_male[,3] <- apply(Gini_male_USA, 1,mean)
Summary_USA_male[,4] <- apply(Gini_male_USA, 1,max)
rownames(Summary_USA_male )<-1959:2020
colnames(Summary_USA_male )<-c("min","median","mean","max")

savepdf("Fig_2b")
matplot(year_USA,Summary_USA_male,type = 'l',lty = c(1,1,1,1),col=c(2,3,4,5),main="Male"
        ,ylab="",xlab="Year",ylim = c(min(Summary_USA_female,Summary_USA_male),max(Summary_USA_female,Summary_USA_male)+0.09))
dev.off()


