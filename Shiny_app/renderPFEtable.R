library(tidyverse)
pfe.data <- readRDS('./pfe_data.Rds')

gen_pfe_table <- function(.country, .division, .gender, .pca) {
  df <- pfe.data %>% 
    filter(country == .country, pol_division == .division, gender == tolower(.gender),
           pca == .pca) %>%
    select(c(-(gender:pol_division), -pca))
  l.table <- df %>% 
    filter(method == "FMP-ANOVA") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select(h, KLD, JSD)
  
  r.table <- df %>% 
    filter(method != "FMP-ANOVA") %>% 
    select(-method) %>% 
    pivot_wider(names_from = metric, values_from = value) %>% 
    select(KLD, JSD)
  
  temp <- cbind(l.table, r.table)
  last.row <- format(round(colMeans(temp), 2), nsmall = 2)
  temp <- temp  %>%
    apply(2, function(x) format(round(x, 2), nsmall = 2)) %>% 
    as.data.frame() 
  temp$h <- as.numeric(temp$h)
  last.row[1] <- "Mean"
  rbind(temp[1:5, ], last.row, temp[-(1:5), ], last.row)
}
