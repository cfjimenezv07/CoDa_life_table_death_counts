library(tidyverse)
pfe.data <- readRDS('./pfe_data.Rds')

gen_pfe_table <- function(.country, .division, .gender, .pca,.fh) {
  df <- pfe.data %>% 
    filter(country == .country, pol_division == .division, gender == tolower(.gender),
           pca == .pca,h== .fh) %>%
    select(c(-(gender:pol_division), -pca,-h))
  
  FMP.table <- df %>% 
    filter(method == "FMP") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select(KLD, JSD)
  
  FM.table <- df %>% 
    filter(method == "FM") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select(KLD, JSD)
  
  TNH.table <- df %>% 
    filter(method == "TNH") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select( KLD, JSD)
  
  GSY.table <- df %>% 
    filter(method == "GSY") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select( KLD, JSD)
  
  MEM.table <- df %>% 
    filter(method == "MEM") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select( KLD, JSD)
  
  temp <- cbind(FMP.table, FM.table,TNH.table, GSY.table,MEM.table)

}

#render the curves for the plots

pfe.curves <- readRDS('./pfe_curves.Rds')

gen_pfe_curves_table <- function(.country, .division, .gender, .pca,.fh) {
  df <- pfe.curves %>% 
    filter(country == .country, pol_division == .division, gender == tolower(.gender),
           pca == .pca,h== .fh) %>%
    select(-gender,-country,-pol_division, -pca,-h)
  
  Actual.table <- df %>% 
    filter(method == "Actual") %>% 
    select(-method)
  
  FMP.table <- df %>% 
    filter(method == "FMP") %>% 
    select(-method) 
  
  
  FM.table <- df %>% 
    filter(method == "FM") %>% 
    select(-method) 
  
  TNH.table <- df %>% 
    filter(method == "TNH") 
  
  GSY.table <- df %>% 
    filter(method == "GSY") 
  
  MEM.table <- df %>% 
    filter(method == "MEM") 
  
  temp <- cbind(Actual.table[,1],FMP.table[,1], FM.table[,1],TNH.table[,1], GSY.table[,1],MEM.table[,1])
  
}
