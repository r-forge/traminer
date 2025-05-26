
unirarcat <- function(bootout, clustering, clusnb, assoc, transformation=FALSE) {
  
  # Ensure correct bootout object
  if(!all(names(bootout) == c("B", "optimal.number", "cluster.solution", 
                              "assoc.char", "original.cluster", "original.assoc", 
                              "coefficients", "errors"))) {
    stop("Please give the output of the function regressboot as argument")
  }
  # Ensure correct size of the clustering solution
  if(nrow(bootout$cluster.solution) != length(clustering)) {
    stop("Please give a clustering solution that has the same size as the original dataset")
  }
  # Ensure that the association can be retrieved
  if(!(assoc %in% bootout$assoc.char)) {
    stop("Please check that the assoc argument has the correct format 
         (as in bootout$assoc.char)")
  }
  # Ensure that the cluster number exists
  if(!(clusnb %in% clustering)) {
    stop("Please give a cluster number as appearing in the clustering solution")
  }
  
  # Retrieve the average marginal effects and standard errors from the bootstrap procedure
  ame <- bootout$coefficients[[assoc]][clustering == clusnb,]
  error <- bootout$errors[[assoc]][clustering == clusnb,]
  
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
  
  # Function from the DescTools package
  if(transformation) prep$ame <- DescTools::FisherZ(prep$ame)
  
  rob <- suppressMessages(lme4::lmer(ame ~ (1|id) + (1|bootstrap), 
                                     weights = prep$stweight, data = prep))
  output <- summary(rob)
  
  if(transformation) {
    
    output$coefficients[1] <- DescTools::FisherZInv(output$coefficients[1])
    output$coefficients[2] <- output$coefficients[2] * (1 - tanh(output$coefficients[1])^2)
    output$varcor$bootstrap <- DescTools::FisherZInv(as.numeric(output$varcor$bootstrap))
    output$varcor$id <- DescTools::FisherZInv(as.numeric(output$varcor$id))
  }
  
  return(list("nobs" = nrow(prep), 
              "pooled.ame" = output$coefficients[1],
              "standard.error" = output$coefficients[2],
              "bootstrap.deviation" = sqrt(as.numeric(output$varcor$bootstrap)),
              "individual.deviation" = sqrt(as.numeric(output$varcor$id)),
              "bootstrap.ranef" = lme4::ranef(rob)$bootstrap[,1],
              "individual.ranef" = lme4::ranef(rob)$id[,1]))
}

