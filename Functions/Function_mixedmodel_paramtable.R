# Function mixedmodel_paramtable: Prepare a table to output a mixed model to a csv table, with added cohen's d effect size.

mixedmodel_paramtable <- function(mymodel) { 
  
  # Extract raw coefficients
  raw_params <- as.data.frame(model_parameters(mymodel))
  
  # Extract standardized coefficients (Cohen's d)
  std_params <- as.data.frame(suppressWarnings(model_parameters(mymodel, standardize = "posthoc", effectsize_type = "cohens_d",  include_random = FALSE)))
  
  # Combine raw and standardized outputs
  paramtable <- left_join(
    raw_params,
    std_params %>% select(Parameter, d = Std_Coefficient),
    by = "Parameter") %>%
    relocate(c("Parameter","Coefficient","SE","z","p","d")) %>%
    mutate(across(c("Coefficient", "SE", "z"), \(x) round(x, 2))) %>%
    mutate(across(c("p"), \(x) round(x, 3))) %>%
    mutate(across(c("d"), \(x) round(x, 3)))
  
  return(paramtable)
    
}
