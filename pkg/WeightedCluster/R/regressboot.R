
regressboot <- function(formula, data, diss, B=500, count=FALSE, 
                        algo="pam", method="ward.D", 
                        fixed=FALSE, kcluster=10, cqi="CH",
                        parallel="no", ncpus=1, cl=NULL) {
  
  # Ensure dissimilarity matrix and dataset have the same size
  if(nrow(data) != nrow(diss)) {
    stop("Please verify that the dissimilarity matrix corresponds to the dataset")
  }
  if(length(kcluster) > 1) {
    stop("Please give a single number as kcluster (see documentation)")
  }
  
  res <- boot::boot(diss, amemat, B, formula=formula, dataset=data, 
                    algo=algo, method=method, cqi=cqi, kcluster=kcluster, fixed=fixed, 
                    count=count, parallel=parallel, ncpus=ncpus, cl=cl)
  originres <- res$t0
  bootres <- t(res$t)
  
  # Create a list of associations
  modelDF <- model.frame(formula, data)
  varlist <- unlist(lapply(colnames(modelDF[-1]), function(i) {
    if(is.factor(modelDF[[i]])) {
      paste0(i, levels(modelDF[[i]])[-1])
    } else {
      if(is.character(modelDF[[i]])) stop("Covariates class should be either factor or numeric")
      i
    }
  }))
  varlist <- sort(varlist)
  n <- nrow(data)
  
  # Extract information once
  first <- bootres[1:n,]
  optimal.number <- apply(first, 2, function(x) length(unique(na.omit(x))))
  cluster.solution <- apply(first, 2, function(x) as.numeric(factor(x), 
                                                             levels = unique(na.omit(x))))
  
  original.assoc <- list()
  effects <- list()
  errors <- list()
  for (i in seq_along(varlist)) {
    
    if(i == 1) {
      effects[[varlist[i]]] <- first
      original.cluster <- as.integer(factor(originres[((i-1)*n+1):(i*n)]))
    } else {
      # The way the function amemat is constructed, the estimates are stacked on top of each other
      # Starting with the AMEs for each association, then the corresponding standard errors
      effects[[varlist[i]]] <- bootres[((i-1)*n+1):(i*n),]
    }
    errors[[varlist[i]]] <- bootres[(length(varlist)*n+1):((length(varlist)+i)*n),]
    original.assoc[[varlist[i]]] <- unique(originres[((i-1)*n+1):(i*n)])
  }
  
  return(list("B" = B, 
              "optimal.kcluster" = optimal.number, 
              "cluster.solution" = cluster.solution, 
              "covar.name" = varlist,
              "original.cluster" = original.cluster,
              "original.ame" = original.assoc,
              "bootstrap.ame" = effects, 
              "std.err" = errors))
}
