
amemat <- function(diss, indices, formula, data1, kmedoid = FALSE, hclust.method="ward.D",
                   fixed = FALSE, ncluster=10, cqi="HC", debug=FALSE) {
  data <- data1
  if(!(cqi %in% c("PBC", "HG", "HGSD", "ASW", "ASWw", "CH", "R2", "CHsq", 
                   "R2sq", "HC"))) {
    stop("Incorrect evaluation measure")
  }
  if(ncluster < 2) stop("At least two clusters required")
  
  # Bootstrap dissimilarities
  d <- diss[indices,indices]
  
  
  ## Clustering the boostrapped sample
  if(fixed||ncluster==2) {
	if(kmedoid) {
		clustering <- wcKMedoid(diss=d, k=ncluster, cluster.only=TRUE)
	} else {
		tree <- fastcluster::hclust(as.dist(d), method = hclust.method)
		clustering <- cutree(tree, k=ncluster)
	}
  }else{ ## Find the best number of group
	if(kmedoid) {
		quality <- wcKMedRange(d, kvals=2:ncluster)
	} else {
		tree <- fastcluster::hclust(as.dist(d), method = hclust.method)
		quality <- as.clustrange(tree, diss=d, ncluster=ncluster)
	}
	if(cqi=="HC"){
		best <- which.min(quality$stats[, cqi])
	}else{
		best <- which.max(quality$stats[, cqi])
	}
	clustering <- quality$clustering[, row.names(quality$stats)[best]]
  }
  clustering <- as.character(clustering)
  ## Association study
  # Bootstrap covariates
  bdata <- data[indices,]
  # Recreate ids for safety
  bdata$id <- indices
  # Changing bootstrap clustering
  bdata[,all.vars(formula)[1]] <- clustering
  
  # Preallocate list
  AMElist <- list()
  
  modelDF <- model.frame(formula, bdata)
  modelFormula <- paste("membership ~", paste(colnames(modelDF)[-1], collapse = " + "))
  for(i in unique(clustering)) {
    
    if(debug) print(i)
    
    modelDF$membership <- modelDF[,1] == i
    mod <- glm(modelFormula, modelDF, family = "binomial")
    AMElist[[i]] <- summary(margins::margins(mod))
  }
  
  output <- createAMEmatrix(clustering, AMElist, indices)
  return(output)
#  m <- match(seq_len(nrow(data)), output[, "id"])
  # With all individuals (and missing observations depending on the sample)
  #output <- dplyr::left_join(data.frame(id = ), output, by = "id")
 # return(c(as.matrix(output[m, -1])))
}

createAMEmatrix <- function(clustering, AMElist, indices=seq_len(length(clustering))){
	 ## ndpul Avoid storing too much duplicate info
	output_list <- vector("list", length(unique(clustering)))

	ndupl <-  !duplicated(indices)
	clustering <- as.numeric(factor(clustering))
	for(i in unique(clustering)){
		# One estimate for each association and each individual in this bootstrap cluster
		ids <- indices[clustering == i & ndupl]
		tmp <- AMElist[[i]]
		effects <- matrix(tmp$AME, nrow = length(ids), ncol = nrow(tmp), byrow = T)
		colnames(effects) <- tmp$factor
		errors <- matrix(tmp$SE, nrow = length(ids), ncol = nrow(tmp), byrow = T)
		colnames(errors) <- paste(tmp$factor, "SE")
		
		output_list[[i]] <- cbind(ids, clustering=i, effects, errors)
	}
	output <- do.call(rbind, output_list)
	m <- match(seq_len(length(clustering)), output[, 1])
	return(as.matrix(output[m, -1]))
}
