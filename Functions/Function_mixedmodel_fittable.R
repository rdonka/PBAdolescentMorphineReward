# Function mixedmodel_fittable: Prepare a table to output the fit stats of a mixed model to a csv table.
# NOTE: This function uses the performance package.

mixedmodel_fittable <- function(fullmodel, nullmodel) { 
  
  modelcomp <- as.data.frame(anova(nullmodel, fullmodel, test = "LRT"))
  
  # Extract metrics
  full_df <- df.residual(fullmodel)
  full_icc <- icc(fullmodel)
  full_r2 <- r2(fullmodel)
  X2 <- modelcomp$`Chisq`[2]
  full_loglik <- logLik(fullmodel)
  full_aic <- AIC(fullmodel)
  full_bic <- BIC(fullmodel)
  
  null_loglik <- logLik(nullmodel)
  null_aic <- AIC(nullmodel)
  null_bic <- BIC(nullmodel)
  
  comp_p <- modelcomp$`Pr(>Chisq)`[2]
  
  # Make a table
  fittable <- data.frame(
    Metric = c("df", "ICC", "Marginal R2", "X2", 
               "Delta AIC", "Delta BIC", "Delta Log-Likelihood", 
               "Model Comparison p",
               "Conditional R2", 
               "Model AIC", "Model BIC", "Model Log-Likelihood", 
               "Null AIC", "Null BIC", "Null Log-Likelihood"),
    Value = c(
      full_df, # Residual degrees of freedom
      full_icc$ICC_adjusted,      # Adjusted ICC
      full_r2$R2_marginal,  # Marginal R2
      X2, # Chi squared
      
      full_aic-null_aic, # Delta AIC
      full_bic-null_bic, # Delta BIC
      full_loglik-null_loglik, # Delta Log-Likelihood
      
      comp_p, # Model comparison p-value
      
      full_r2$R2_conditional, # Conditional R2
      
      full_aic, # Full model AIC
      full_bic, # Full model BIC
      full_loglik, # Full model Log-Likelihood
      
      null_aic, # Null model AIC
      null_bic, # Null model BIC
      null_loglik # Null model Log-Likelihood
    )
  )

  fittable$Value <- round(fittable$Value, 3)
  
  return(fittable)
  
  }
