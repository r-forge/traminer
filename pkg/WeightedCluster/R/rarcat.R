
rarcat <- function(diss, covar, df, 
                   clustering=NULL, robust=TRUE, B=500, count=FALSE, 
                   algo="pam", method="ward.D", 
                   fixed=FALSE, ncluster=10, eval="CH",
                   parallel="no", ncpus=1, cl=NULL,
                   transformation=FALSE, conflevel=0.05, digits=3) {
  
  # Ensure that the associations can be retrieved
  if(!all(covar %in% colnames(df))) {
    stop("Please check that the covariates names correspond to columns in the dataset")
  }
  
  if(!robust & is.null(clustering)) {
    stop("Bootstrap procedure not asked so please input a clustering solution")
  }
  
  # Run the bootstrap procedure if asked
  bootout <- NULL
  if(robust) {
    bootout <- regressboot(diss, covar, df, B=B, count=count, 
                           algo=algo, method=method, 
                           fixed=fixed, ncluster=ncluster, eval=eval,
                           parallel=parallel, ncpus=ncpus, cl=cl)
    clustering <- bootout$original.cluster
  }
  
  # Ensure correct size of the clustering solution
  if(nrow(df) != length(clustering)) {
    stop("Please give a clustering solution that has the same size as the original dataset")
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
        res <- unirarcat(bootout, clustering, i, robust$factor[j], transformation)
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
