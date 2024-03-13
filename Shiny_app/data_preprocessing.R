#Making data ready for the pfe table

#gettig the names regarding the departments or states
names.france <- sapply(readRDS("names/names_departments.rds"), function(x) x)
names.us <- sapply(readRDS("names/names_states.rds"), function(x) x)
names <- list(France = names.france, USA = names.us)

get_pfe_dataset <- function(.country, .method, .metric, .gender, .pca) {
  #reading the data
  path <- paste0("./../datasets_shiny_app copy/PFE/",.country, "/", 
                 paste(.method, .metric, .gender, .pca, sep = "_"), ".rds")
  df <- readRDS(path)
  colnames(df) <- names[[.country]] #labeling the departments or states
  df <- as.data.frame(df)
  #adding the identifiers such
  df$method <- .method
  df$pca <- .pca
  df$metric <- .metric
  df$h <- 1:10
  df$gender <- .gender
  df$country <- .country
  #pivoting dataset
  df %>% 
    pivot_longer(-(method:country), names_to = "pol_division")
}

#generating data
country.list <- as.list(rep(c("France", "US"), each = 8))
method.list <- as.list(rep(c("FM", "FMP"), each = 8))
metric.list <- as.list(rep(c("JSD", "KLD"), each = 8))
gender.list <- as.list(rep(c("female", "male"), each = 8))
pca.list <- as.list(rep(c("EVR", "K"), each = 8))


lista <- expand_grid(country = c("USA"), method = c("FM", "FMP"),
            metric = c("JSD", "KLD"), gender = c("male", "female"), pca = c("EVR", "K")) %>%
apply(2, function(x) as.list(x))

pfe.result <- mapply(get_pfe_dataset, .country = lista$country, .method = lista$method, 
       .metric = lista$metric, .gender = lista$gender, .pca = lista$pca, SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .) %>%
  mutate(value = round(value*100, 2), method = paste0(method, "-ANOVA"))

View(pfe.result)
saveRDS(pfe.result, "./pfe_data.rds")

#ife dataset

aux <- function(i, df, .country, .gender, .method, .coverage, .pca) {
  temp <- df[i, , ]
  temp[, -3] <- round(temp[, -3]*100, 2)
  temp[, 3] <- round(temp[, 3]*1e-3, 2)
  colnames(temp) <- c("ECP", "CPD", "IS")
  temp <- as.data.frame(temp)
  temp$method <- .method 
  temp$pca <- .pca
  temp$h <- 1:10
  temp$coverage <- .coverage
  temp$gender <- .gender
  temp$country <- .country
  temp$pol_division = names[[.country]][i]
  temp
}

get_ife_dataset <- function(.country, .method, .gender, .coverage, .pca) {
  #reading the data
  path <- paste0("./../datasets_shiny_app copy/IFE/",.country,"/", 
                 paste("Emp_cov", .gender, .coverage, sep = "_"), 
                 ifelse(.method == "FM", "_means_", "_"),
                 paste(.pca, "1lrc", sep = "_"), ".rds")
  df <- readRDS(path)
  do.call(rbind.data.frame, 
               lapply(1:51, aux, df, .country, .gender, .method, .coverage, .pca))
}

lista <- expand_grid(country = c("USA"), method = c("FM", "FMP"),
                     gender = c("male", "female"), coverage = c(80, 95), pca = c("EVR", "K")) %>%
  apply(2, function(x) as.list(x))

ife.result <- mapply(get_ife_dataset, .country = lista$country, .method = lista$method, 
       .gender = lista$gender, .coverage = lista$coverage, 
       .pca = lista$pca, SIMPLIFY = F) %>%
  do.call(rbind.data.frame, .) %>%
  mutate(method = paste0(method, "-ANOVA")) %>%
  pivot_longer(ECP:IS, names_to = "metric") %>%
  relocate(metric, everything())
saveRDS(ife.result, "./ife_data.rds")

View(ife.result)
