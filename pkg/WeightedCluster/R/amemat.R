
amemat <- function(diss, indices, df, covar, algo="pam", method="ward.D",
                   fixed = FALSE, ncluster=10, eval="CH", count=FALSE) {
  
  if(!(eval %in% c("PBC", "HG", "HGSD", "ASW", "ASWw", "CH", "R2", "CHsq", 
                   "R2sq", "HC"))) {
    stop("Incorrect evaluation measure")
  }
  if(ncluster < 2) stop("At least two clusters required")
  
  # Bootstrap dissimilarities
  d <- diss[indices,indices]
  
  if(algo == "pam") {
    
    k <- 2:ncluster
    quality <- wcKMedRange(d, kvals=k)
    
  } else if(algo == "hierarchical") {
    
    tree <- fastcluster::hclust(as.dist(d), method = method)
    if(ncluster == 2) {
      clustering <- cutree(tree, k=ncluster)
    } else {
      quality <- WeightedCluster::as.clustrange(tree, diss=d, ncluster=ncluster)
    }
  } else {
    
    stop("Clustering algorithm currently not supported")
  }
  
  if(!fixed) {
    quality$stats <- arrange(quality$stats, !!as.symbol(eval))
    if(eval == "HC") quality$stats <- arrange(quality$stats, desc(!!as.symbol(eval)))
  }
  
  if(ncluster > 2 | algo == "pam") {
    clustering <- quality$clustering[,row.names(quality$stats)[nrow(quality$stats)]]
  }
  
  # Bootstrap covariates
  replicate <- df[indices,]
  # Recreate ids for safety
  replicate$id <- indices
  # Optimal solution
  replicate$solution <- clustering
  
  # Create formula object
  formula <- paste("membership ~", paste(covar, collapse = " + "))
  # Preallocate list
  output_list <- vector("list", length(unique(replicate$solution)))
  
  for(i in unique(replicate$solution)) {
    
    if(count) print(i)
    
    replicate$membership <- replicate$solution == i
    mod <- glm(formula, replicate, family = "binomial")
    tmp <- summary(margins::margins(mod))
    
    # One estimate for each association and each individual in this bootstrap cluster
    sub <- replicate %>% 
      filter(solution == i, !duplicated(id)) %>% 
      select(id)
    effects <- matrix(tmp$AME, nrow = nrow(sub), ncol = nrow(tmp), byrow = T)
    colnames(effects) <- tmp$factor
    errors <- matrix(tmp$SE, nrow = nrow(sub), ncol = nrow(tmp), byrow = T)
    colnames(errors) <- paste(tmp$factor, "SE")
    
    output_list[[i]] <- cbind(sub, effects, errors)
  }
  
  output <- do.call(rbind, output_list)
  # With all individuals (and missing observations depending on the sample)
  output <- left_join(data.frame(id = seq_len(nrow(df))), output, by = "id")
  return(c(as.matrix(output[,-1])))
}
