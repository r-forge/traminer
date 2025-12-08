
bootpool <- function(bootout, clustering, clusnb, covar, fisher_transform = FALSE) {
  
  # Ensure correct bootout object
  if(!all(names(bootout) == c("B", "optimal.kcluster", "cluster.solution", 
                              "covar.name", "original.cluster", "original.ame", 
                              "bootstrap.ame", "std.err"))) {
    stop(" [!] Please give the output of the function regressboot as argument")
  }
  # Ensure correct size of the clustering solution
  if(nrow(bootout$cluster.solution) != length(clustering)) {
    stop(" [!] Please give a clustering solution that has the same size as the original dataset")
  }
  # Ensure that the association can be retrieved
  if(!(covar %in% bootout$covar.name)) {
    stop(" [!] Please check that the covar argument has the correct format 
         (as in bootout$covar.name)")
  }
  # Ensure that the cluster number exists
  if(!(clusnb %in% clustering)) {
    stop(" [!] Please give a cluster number as appearing in the clustering solution")
  }
  
  # Retrieve the average marginal effects and standard errors from the bootstrap procedure
  ame <- bootout$bootstrap.ame[[covar]][clustering == clusnb,]
  error <- bootout$std.err[[covar]][clustering == clusnb,]
  
  # Prepare the data to input in the multilevel model
  tmp <- data.frame(clustering = clustering, id = seq(1, length(clustering)))
  prep <- data.frame(bootstrap = rep(1:ncol(ame), nrow(ame)),
                     id = rep(dplyr::filter(tmp, clustering == clusnb)$id, each=ncol(ame)),
                     ame = c(t(ame)), 
                     sterror = c(t(error)))
  prep <- dplyr::filter(prep, !is.na(ame))
  # Weights based on the standard errors
  prep$weight <- 1/prep$sterror^2
  prep$stweight <- prep$weight/mean(prep$weight)
  
  # Atanh for Fisher's z-transformation
  if(fisher_transform) prep$ame <- atanh(prep$ame)
  
  rob <- suppressMessages(lme4::lmer(ame ~ (1|id) + (1|bootstrap), 
                                     weights = prep$stweight, data = prep))
  output <- summary(rob)
  
  if(fisher_transform) {
    #tanh for inverse Fisher's z-transformation
    output$coefficients[1] <- tanh(output$coefficients[1])
    output$coefficients[2] <- output$coefficients[2] * (1 - tanh(output$coefficients[1])^2)
    output$varcor$bootstrap <- tanh(as.numeric(output$varcor$bootstrap))
    output$varcor$id <- tanh(as.numeric(output$varcor$id))
  }
  
  return(list("nobs" = nrow(prep), 
              "pooled.ame" = output$coefficients[1],
              "standard.error" = output$coefficients[2],
              "bootstrap.stddev" = sqrt(as.numeric(output$varcor$bootstrap)),
              "individual.stddev" = sqrt(as.numeric(output$varcor$id)),
              "bootstrap.ranef" = lme4::ranef(rob)$bootstrap[,1],
              "individual.ranef" = lme4::ranef(rob)$id[,1]))
}
