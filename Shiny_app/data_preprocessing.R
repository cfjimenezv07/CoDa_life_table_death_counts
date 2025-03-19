#Making data ready for the pfe table

#gettig the names regarding the departments or states
names.france <- sapply(readRDS("names/names_departments.rds"), function(x) x)
names.us <- sapply(readRDS("names/names_states.rds"), function(x) x)
names <- list(USA = names.us)

get_pfe_dataset <- function(.country, .method, .metric, .gender, .pca) {
  #reading the data
  path <- paste0("./datasets_shiny_app/PFE/",.country, "/", 
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
country.list <- as.list(rep(c("France", "USA"), each = 8))
method.list <- as.list(rep(c("FM", "FMP","TNH","GSY","MEM"), each = 8))
metric.list <- as.list(rep(c("JSD", "KLD"), each = 8))
gender.list <- as.list(rep(c("female", "male"), each = 8))
pca.list <- as.list(rep(c("EVR", "K"), each = 8))


lista <- expand_grid(country = c("USA"), method = c("FM", "FMP","TNH","GSY","MEM"),
            metric = c("JSD", "KLD"), gender = c("male", "female"), pca = c("EVR", "K")) %>%
apply(2, function(x) as.list(x))

pfe.result <- mapply(get_pfe_dataset, .country = lista$country, .method = lista$method, 
       .metric = lista$metric, .gender = lista$gender, .pca = lista$pca, SIMPLIFY = F) %>% 
  do.call(rbind.data.frame, .) %>%
  mutate(value = round(value*100, 2), method = method)

View(pfe.result)
saveRDS(pfe.result, "./pfe_data.rds")



##############################
# Datasets to get the actual forecasts
##############################


get_pfe_curves <- function(.method, .gender, .pca) {
  # Define the path to the dataset (without changing the path)
  path <- paste0("./datasets_shiny_app/PFE/USA/Results_curves/", 
                 paste(.method, .gender, .pca, sep = "_"), ".rds")
  
  # Read the dataset
  df <- readRDS(path)
  
  # Define the state names (since the dataset has 51 states)
  state_names <- c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", 
                   "Delaware", "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", 
                   "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota", 
                   "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey", 
                   "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", 
                   "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", 
                   "Vermont", "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming")
  
  # Initialize an empty list to store results
  result_list <- list()
  
  # Iterate over each state (i.e., the first dimension of the array)
  for (i in 1:51) {
    # Extract data for the current state (state `i`)
    state_data <- df[[i]]
    
    # Extract values for the given forecast horizon (h)
    # Let's assume the column for `h` is the first column for h=1, second for h=2, etc.
    # You can adjust this depending on your actual data structure
    for (h_value in 1:10) {  # Assuming h ranges from 1 to 10
      df_state <- state_data[, h_value]
      
      # Convert to a data frame and add the necessary metadata
      df_long <- data.frame(
        value = df_state,
        method = .method,
        gender = .gender,
        pca = .pca,
        pol_division = state_names[i],
        h = h_value
      )
      
      # Add this data frame to the result list
      result_list[[length(result_list) + 1]] <- df_long
    }
  }
  
  # Combine all the data into a single data frame
  final_df <- do.call(rbind, result_list)
  
  return(final_df)
}

# Example usage with multiple combinations of method, gender, and pca
lista2 <- expand_grid(method = c("Actual","FM", "FMP", "TNH", "GSY", "MEM"),
                      gender = c("male", "female"),
                      pca = c("EVR", "K"))

# Get all the data by applying the function to all combinations
pfe.curves <- mapply(get_pfe_curves, 
                     .method = lista2$method,
                     .gender = lista2$gender,
                     .pca = lista2$pca, 
                     SIMPLIFY = FALSE) %>% 
  do.call(rbind.data.frame, .) %>%
  mutate(value = round(value * 100, 2))  # If needed, scale the values
pfe.curves$country <- "USA"
# View the resulting dataframe
View(pfe.curves)
saveRDS(pfe.curves, "./pfe_curves.rds")


