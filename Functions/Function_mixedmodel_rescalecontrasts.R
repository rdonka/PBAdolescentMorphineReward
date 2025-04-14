# Function mixedmodel_rescalecontrasts: Rescale the coefficients, SEs, and CIs in a model table to raw units.
# Used with models fit to scaled data. Multiplies scaled values by the raw unit SD and adjusts the intercept
# so that interpretation of the model outputs is in raw units.

mixedmodel_rescalecontrasts <- function(mycontrasts, myrawsd) { 
  
  mycontrasts$estimate <- round((mycontrasts$estimate * myrawsd),2) # Adjust coefficients for scaling
  mycontrasts$SE <- round((mycontrasts$SE * myrawsd),2) # Adjust SEs for scaling

  return(mycontrasts)
}
