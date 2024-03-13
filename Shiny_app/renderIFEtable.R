library(tidyverse)
ife.data <- readRDS("./ife_data.rds")

gen_ife_table <- function(.country, .divison, .gender, .coverage, .pca) {
  df <- ife.data %>% 
    filter(country == .country, pol_division == .divison, 
           gender == tolower(.gender), coverage == as.numeric(gsub("%", "", .coverage)),
           pca == .pca) %>%
    select(c(-(gender:pol_division), -pca))
  l.table <- df %>% 
    filter(method == "FMP-ANOVA") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select(h, ECP, CPD, IS) %>% 
    round(2)
  
  r.table <- df %>% 
    filter(method != "FMP-ANOVA") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select(ECP, CPD, IS) %>% 
    round(2)
  
  temp <- cbind(l.table, r.table)
  last.row <- format(round(colMeans(temp), 2), nsmall = 2)
  temp <- temp  %>%
    apply(2, function(x) format(round(x, 2), nsmall = 2)) %>% 
    as.data.frame()
  temp$h <- as.numeric(temp$h)
  last.row[1] <- "Mean"
  rbind(temp[1:5, ], last.row, temp[-(1:5), ], last.row)
}
