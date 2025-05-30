
rarcat <- function(formula, data, diss, 
                   robust=TRUE, B=500, count=FALSE, 
                   algo="pam", method="ward.D", 
                   fixed=FALSE, kcluster=10, cqi="CH",
                   parallel="no", ncpus=1, cl=NULL,
                   fisher_transform=FALSE, conflevel=0.05, digits=3) {
  
  # Run the bootstrap procedure if asked
  bootout <- NULL
  if(robust) {
    bootout <- regressboot(formula, data, diss, B=B, count=count, 
                           algo=algo, method=method, 
                           fixed=fixed, kcluster=kcluster, cqi=cqi,
                           parallel=parallel, ncpus=ncpus, cl=cl)
  }
  
  modelDF <- model.frame(formula, data)
  modelFormula <- paste("membership ~", paste(colnames(modelDF)[-1], collapse = " + "))
  
  original <- data.frame()
  assess <- data.frame()
  for(i in unique(modelDF[,1])) {
    
    modelDF$membership <- modelDF[,1] == i
    mod <- glm(modelFormula, modelDF, family = "binomial")
    tmp <- summary(margins::margins(mod))
    #print(tmp)
    
    # Set up with the specific associations
    if(ncol(original) == 0) {
      original <- rbind(original, data.frame(factor = tmp$factor))
      assess <- rbind(assess, data.frame(factor = tmp$factor))
    } 
    
    # Average marginal effects with confidence intervals
    original[,paste("cluster", i)] <- round(tmp$AME, digits)
    original[,paste0("c", i, " lower")] <- 
      round(tmp$AME - qnorm(1 - conflevel/2)*tmp$SE, digits)
    original[,paste0("c", i, " upper")] <- 
      round(upper <- tmp$AME + qnorm(1 - conflevel/2)*tmp$SE, digits)
    
    if(!is.null(bootout)) {
      
      assess[,paste("cluster", i)] <- NA
      assess[,paste0("c", i, " lower")] <- NA
      assess[,paste0("c", i, " upper")] <- NA
      for(j in 1:nrow(assess)) {
        
        # Run the pooling function for each combination of cluster and covariate
        res <- bootpool(bootout, modelDF[,1], i, assess$factor[j], fisher_transform)
        # Formula for the prediction interval
        var <- qt(1 - conflevel/2, bootout$B - 2)*
          sqrt(res$standard.error^2 + res$bootstrap.stddev^2)
        
        assess[j,paste("cluster", i)] <- round(res$pooled.ame, digits)
        assess[j,paste0("c", i, " lower")] <- round(res$pooled.ame - var, digits)
        assess[j,paste0("c", i, " upper")] <- round(res$pooled.ame + var, digits)
      }
    }
  }
  
  return(list("original.analysis" = original,
              "robust.analysis" = assess))
}
