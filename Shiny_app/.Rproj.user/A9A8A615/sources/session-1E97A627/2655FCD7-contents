# Tables of comparison between lrc and cov with ARIMA.
setwd("~/My Drive/Fall 2023/STAT 397/PhD project 4/Revisions from JCGS/Rcodes/")

# USA
All_errors_tables_USA_static <- readRDS("./All_errors_tables_USA_static.rds")
All_errors_tables_USA_FM_static <- readRDS("./All_errors_tables_USA_FM_static.rds")

All_errors_table_USA_dynamic <- readRDS("./All_errors_table_USA_dynamic.rds")
All_errors_table_USA_FM_dynamic <- readRDS("./All_errors_table_USA_FM_dynamic.rds")



####################################################################################################
# Tables
####################################################################################################

#1. USA
USA_results_FMP <- matrix(NA,nrow=10,ncol = 8)
USA_results_FM <- matrix(NA,nrow=10,ncol = 8)
for (i in 1:4) {
  USA_results_FMP[,(2*i-1)]=apply(All_errors_table_USA_dynamic[[i]], 1, mean)
  USA_results_FMP[,2*i]=apply(All_errors_tables_USA_static[[i]], 1, mean)
  USA_results_FM[,(2*i-1)]=apply(All_errors_table_USA_FM_dynamic[[i]], 1, mean)
  USA_results_FM[,2*i]=apply(All_errors_tables_USA_FM_static[[i]], 1, mean)
}
library(xtable)
xtable(USA_results_FMP)
xtable(USA_results_FM)

mean_USA_FMP <- round(apply(USA_results_FMP , 2, mean),2)
mean_USA_FM <- round(apply(USA_results_FM , 2, mean),2)
####################################################################################################
# France
All_errors_tables_France_static <- readRDS("./All_erros_table_France_static.rds")
All_errors_tables_France_FM_static <- readRDS("./All_errors_table_France_FM_static.rds")

All_errors_table_France_dynamic <- readRDS("./All_errors_table_France_dynamic.rds")
All_errors_table_France_FM_dynamic <- readRDS("./All_errors_table_France_FM_dynamic.rds")

#2. France
France_results_FMP <- matrix(NA,nrow=10,ncol = 8)
France_results_FM <- matrix(NA,nrow=10,ncol = 8)
for (i in 1:4) {
  France_results_FMP[,(2*i-1)]=apply(All_errors_table_France_dynamic[[i]], 1, mean)
  France_results_FMP[,2*i]=apply(All_errors_tables_France_static[[i]], 1, mean)
  France_results_FM[,(2*i-1)]=apply(All_errors_table_France_FM_dynamic[[i]], 1, mean)
  France_results_FM[,2*i]=apply(All_errors_tables_France_FM_static[[i]], 1, mean)
}

xtable(France_results_FMP)
xtable(France_results_FM)


mean_France_FMP <- round(apply(France_results_FMP , 2, mean),2)
mean_France_FM <- round(apply(France_results_FM , 2, mean),2)
