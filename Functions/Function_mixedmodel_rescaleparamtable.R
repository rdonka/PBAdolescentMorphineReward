# Function mixedmodel_rescaleparamtable: Rescale the coefficients, SEs, and CIs in a model table to raw units.
# Used with models fit to scaled data. Multiplies scaled values by the raw unit SD and adjusts the intercept
# so that interpretation of the model outputs is in raw units.

mixedmodel_rescaleparamtable <- function(myparamtable, myrawsd, myrawmean) { 
  
  myparamtable$Coefficient <- round((myparamtable$Coefficient * myrawsd),2) # Adjust coefficients for scaling
  myparamtable$Coefficient[which(myparamtable$Parameter == "(Intercept)")] <- round((myparamtable$Coefficient[which(myparamtable$Parameter == "(Intercept)")] + myrawmean),2) # Adjust intercept for scaling
  
  myparamtable$SE <-  round((myparamtable$SE * myrawsd),2) # Adjust coefficients for scaling
  myparamtable$CI_low <-  round((myparamtable$CI_low * myrawsd),2) # Adjust coefficients for scaling
  myparamtable$CI_high <-  round((myparamtable$CI_high * myrawsd),2) # Adjust coefficients for scaling
  
  return(myparamtable)
}
