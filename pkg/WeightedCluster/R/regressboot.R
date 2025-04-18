
regressboot <- function(diss, covar, df, B=500, count=FALSE, 
                        algo="pam", method="ward.D", 
                        fixed=FALSE, ncluster=10, eval="CH",
                        parallel="no", ncpus=1, cl=NULL) {
  
  # Ensure dissimilarity matrix and dataset have the same size
  if(nrow(df) != nrow(diss)) {
    stop("Please verify that the dissimilarity matrix corresponds to the dataset")
  }
  # Ensure that the associations can be retrieved
  if(!all(covar %in% colnames(df))) {
    stop("Please check that the covariates names correspond to columns in the dataset")
  }
  if(length(ncluster) > 1) {
    stop("Please give a single number as ncluster (see documentation)")
  }
  
  res <- boot::boot(diss, amemat, B, df=df, covar=covar, algo=algo, method=method, 
                    eval=eval, ncluster=ncluster, fixed=fixed, count=count, 
                    parallel=parallel, ncpus=ncpus, cl=cl)
  originres <- res$t0
  bootres <- t(res$t)
  
  # Create a list of associations
  varlist <- unlist(lapply(covar, function(i) {
    if(is.factor(df[[i]])) {
      paste0(i, levels(df[[i]])[-1])
    } else {
      i
    }
  }))
  varlist <- sort(varlist)
  
  # Extract information once
  first <- bootres[1:nrow(df),]
  optimal.number <- apply(first, 2, function(x) length(unique(na.omit(x))))
  cluster.solution <- apply(first, 2, function(x) as.numeric(factor(x), 
                                                             levels = unique(na.omit(x))))
  
  original <- list()
  effects <- list()
  errors <- list()
  for (i in seq_along(varlist)) {
    
    if(i == 1) {
      effects[[varlist[i]]] <- first
    } else {
      # The way the function amemat is constructed, the estimates are stacked on top of each other
      # Starting with the AMEs for each association, then the corresponding standard errors
      effects[[varlist[i]]] <- bootres[((i-1)*nrow(df)+1):(i*nrow(df)),]
    }
    errors[[varlist[i]]] <- bootres[(length(varlist)*nrow(df)+1):((length(varlist)+i)*nrow(df)),]
    original[[varlist[i]]] <- unique(originres[((i-1)*nrow(df)+1):(i*nrow(df))])
  }
  
  return(list("B" = B, 
              "optimal.number" = optimal.number, 
              "cluster.solution" = cluster.solution, 
              "assoc.char" = varlist,
              "original.result" = original,
              "coefficients" = effects, 
              "errors" = errors))
}
