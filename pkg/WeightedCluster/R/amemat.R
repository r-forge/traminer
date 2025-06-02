
amemat <- function(diss, indices, formula, dataset, algo="pam", method="ward.D",
                   fixed = FALSE, kcluster=10, cqi="CH", count=FALSE) {
  
  if(!(cqi %in% c("PBC", "HG", "HGSD", "ASW", "ASWw", "CH", "R2", "CHsq", 
                   "R2sq", "HC"))) {
    stop("Incorrect evaluation measure")
  }
  if(kcluster < 2) stop("At least two clusters required")
  
  # Bootstrap dissimilarities
  d <- diss[indices,indices]
  
  if(algo == "pam") {
    
    k <- 2:kcluster
    quality <- wcKMedRange(d, kvals=k)
    
  } else if(algo == "hierarchical") {
    
    tree <- fastcluster::hclust(as.dist(d), method = method)
    if(kcluster == 2) {
      clustering <- cutree(tree, k=kcluster)
    } else {
      quality <- WeightedCluster::as.clustrange(tree, diss=d, ncluster=kcluster)
    }
  } else {
    
    stop("Clustering algorithm currently not supported")
  }
  
  if(!fixed) {
    quality$stats <- dplyr::arrange(quality$stats, !!as.symbol(cqi))
    if(cqi == "HC") quality$stats <- dplyr::arrange(quality$stats, desc(!!as.symbol(cqi)))
  }
  
  if(kcluster > 2 | algo == "pam") {
    clustering <- quality$clustering[,row.names(quality$stats)[nrow(quality$stats)]]
  }
  
  # Bootstrap covariates
  replicate <- dataset[indices,]
  # Recreate ids for safety
  replicate$id <- indices
  # Optimal solution
  replicate[,all.vars(formula)[1]] <- clustering
  
  # Preallocate list
  output_list <- vector("list", length(unique(clustering)))
  
  modelDF <- model.frame(formula, replicate)
  modelFormula <- paste("membership ~", paste(colnames(modelDF)[-1], collapse = " + "))
  
  for(i in unique(clustering)) {
    
    if(count) print(i)
    
    modelDF$membership <- modelDF[,1] == i
    mod <- glm(modelFormula, modelDF, family = "binomial")
    tmp <- summary(margins::margins(mod))
    
    # One estimate for each association and each individual in this bootstrap cluster
    sub <- replicate[replicate[,all.vars(formula)[1]] == i & 
                       !duplicated(replicate$id), "id", drop = FALSE]
    effects <- matrix(tmp$AME, nrow = nrow(sub), ncol = nrow(tmp), byrow = T)
    colnames(effects) <- tmp$factor
    errors <- matrix(tmp$SE, nrow = nrow(sub), ncol = nrow(tmp), byrow = T)
    colnames(errors) <- paste(tmp$factor, "SE")
    
    output_list[[i]] <- cbind(sub, effects, errors)
  }
  
  output <- do.call(rbind, output_list)
  # With all individuals (and missing observations depending on the sample)
  output <- dplyr::left_join(data.frame(id = seq_len(nrow(dataset))), output, by = "id")
  return(c(as.matrix(output[,-1])))
}
