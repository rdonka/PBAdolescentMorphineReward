# Function mixedmodel_paramtable: Prepare a table to output a mixed model to a csv table, with added cohen's d effect size.

mixedmodel_emmeanstable_SessionTime <- function(myemmeans, mymodel) { 
  
  ES <- eff_size(
    myemmeans$emmeans,
    sigma = sigma(mymodel),   # the residual SD from your model
    edf   = Inf                              # or a suitable degrees-of-freedom estimate
  ) 
  
  EStable <- as.data.frame(ES) %>% dplyr::select(contrast, SessionTime, effect.size)
  
  # Extract contrasts# Extract contraeffectsize()sts
  emmeanstable <- as.data.frame(myemmeans$contrasts) %>% 
    left_join(EStable, by=c('contrast','SessionTime')) %>%
    relocate(c("contrast","df","t.ratio","p.value", "effect.size")) %>%
    mutate(across(c("t.ratio"), \(x) round(x, 2))) %>%
    mutate(across(c("p.value"), \(x) round(x, 3))) %>%
    mutate(across(c("effect.size"), \(x) round(x, 3)))
  
  return(emmeanstable)
    
}
