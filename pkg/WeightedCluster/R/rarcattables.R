
rarcattables <- function(clustering, covar, df, bootout=NULL, transformation=FALSE, 
                         conflevel=0.05, digits=3) {
  
  # Ensure correct size of the clustering solution
  if(nrow(df) != length(clustering)) {
    stop("Please give a clustering solution that has the same size as the original dataset")
  }
  # Ensure that the associations can be retrieved
  if(!all(covar %in% colnames(df))) {
    stop("Please check that the covariates names correspond to columns in the dataset")
  }
  
  formula <- paste("membership ~", paste(covar, collapse = " + "))
  
  original <- data.frame()
  robust <- data.frame()
  for(i in unique(clustering)) {
    
    df$membership <- clustering == i
    mod <- glm(formula, df, family = "binomial")
    tmp <- summary(margins::margins(mod))
    #print(tmp)
    
    # Set up with the specific associations
    if(ncol(original) == 0) {
      original <- rbind(original, data.frame(factor = tmp$factor))
      robust <- rbind(robust, data.frame(factor = tmp$factor))
    } 
    
    # Average marginal effects with confidence intervals
    lower <- tmp$AME - qnorm(1 - conflevel/2)*tmp$SE
    upper <- tmp$AME + qnorm(1 - conflevel/2)*tmp$SE
    original[,paste("cluster", i)] <- paste0(round(tmp$AME, digits), " [", 
                                             round(lower, digits), ", ", 
                                             round(upper, digits), "]")
    
    if(!is.null(bootout)) {
      
      robust[,paste("cluster", i)] <- NA
      for(j in 1:nrow(robust)) {
        
        # Run the RARCAT function for each combination of cluster and covariate
        res <- rarcat(bootout, clustering, i, robust$factor[j], transformation)
        # Formula for the prediction interval
        var <- qt(1 - conflevel/2, bootout$B - 2)*
          sqrt(res$standard.error^2 + res$bootstrap.deviation^2)
        min <- res$pooled.ame - var
        max <- res$pooled.ame + var
        
        robust[j,paste("cluster", i)] <- paste0(round(res$pooled.ame, digits), " [", 
                                                round(min, digits), ", ", 
                                                round(max, digits), "]") 
      }
    }
  }
  
  return(list("original.analysis" = original,
              "robust.analysis" = robust))
}
