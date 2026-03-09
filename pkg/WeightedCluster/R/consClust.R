
## clueboot is mostly internal
cl_class_ids.cluebootkmed <- function(x) as.cl_class_ids(x) # functions needed by clueInterface to run 

is.cl_partition.cluebootkmed <- function(x) TRUE

is.cl_hard_partition.cluebootkmed <- function(x) TRUE

clueboot <- function(diss,
                     k, 
                     R=100,
                     int.weights = TRUE,
                     method="dirichlet", 
                     base.clust = "pam"){
  cluster <- base.clust
  n <- nrow(diss)
  ret <- list()
  if(method=="dirichlet"){
    allweights <- matrix(rexp(n * R, 1) , ncol = n, byrow = TRUE)
  }else if(method=="uniform"){
    allweights <- matrix(runif(n * R) , ncol = n, byrow = TRUE)
  }
  
  allweights[allweights<=0] <- min(allweights[allweights>0])/2 ##Avoid 0 weights
  allweights <- n*allweights / rowSums(allweights)
  if(int.weights){
    allweights <- allweights/min(allweights) - 1
    allweights <- round(allweights*(2e+09 - 2)/max(allweights))+1
    allweights <- as.integer(allweights)
    dim(allweights) <- c(R, n)
  }
  
  ## Ensure we have t0
  allweights[1, ] <- 1
  
  diss <- as.dist(diss)
  clustAlgo <- rep(base.clust, R/length(base.clust))
  if(length(clustAlgo) < R){
    dif <- R - length(clustAlgo)
    
    for(i in 1:dif){
      clustAlgo[(R-dif) + i] <- base.clust[i]   
    }
  }
  else{
    
  }
  for(i in 1:R){	## Consider parallel?
    if(clustAlgo[i] %in% c("pam", "single", "complete", "average", "mcquitty",
                             "ward.D", "ward.D2", "centroid" , "median") ){ 
      if(clustAlgo[i] == "pam" ){
        ret[[i]] <- wcKMedoids(diss, k=k, weights=allweights[i, ], cluster.only=TRUE)
      }
      else{
         hc <- fastcluster::hclust(diss, method=clustAlgo[i], members=allweights[i, ])
         ret[[i]] <- cutree(hc, k=k) 
      }
      
    }
    else{
      cat(paste0('[>] ', 
                 clustAlgo[i], 
                 " is not a supported clustering algorithm"),
          sep = "\n")
    }
    
    class(ret[[i]]) <- "cluebootkmed"
  }
  names(ret) <- paste0("part_", 1:R, "_", clustAlgo)
  return(cl_ensemble(list=ret))
}



#-------------------------------------------------------------------------------

consClust <- function(diss, base.clust = "pam", R =100, 
                      kvals = 2:15, cons.method = "SE",  
                      membership = "crisp", k.fixed = TRUE, 
                      agg.method = "cRand", keep.ensemble = TRUE, 
                      parallel = FALSE, progressbar = TRUE){
  tStart <- Sys.time()
  
  cat(paste0("[>] Performing consensus clustering on ", 
             R, " partitions, using: ", paste0(base.clust, collapse = ", ")),
      sep = "\n")
  
  if(parallel){
    message(" [>] Initializing parallel processing...", appendLF = FALSE)
    oplan <- plan(multisession)
    on.exit(plan(oplan), add = TRUE)
    message(" OK.")
  }
  else{
  }
  if (progressbar) {
    if (requireNamespace("progress", quietly = TRUE)) {
      old_handlers <- handlers(handler_progress(format = "(:spin) [:bar] :percent | Elapsed: :elapsed | ETA: :eta | :message"))
      if (!is.null(old_handlers)) {
        on.exit(handlers(old_handlers), add = TRUE)
      }
    }
    else {
      message(" [>] Install the progress package to see estimated remaining time.")
    }
    oldglobal <- handlers(global = TRUE)
    if (!is.null(oldglobal)) {
      on.exit(handlers(global = oldglobal), add = TRUE)
    }
 
  }
  p <- progressor(max(kvals))
  gc()
  it <- as.data.frame(expand.grid(1:R, kvals), col.names = c("R", "j"))
  ff <- foreach(j = kvals, 
                .options.future = list(seed = TRUE, 
                                       packages = c("clue", 
                                                    "fastcluster", 
                                                    "WeightedCluster"),
                                       globals = structure(TRUE, 
                                                add = c("base.clust", "k", "R",
                                              "cons.method")))) %dofuture% {
		#source("cacheFunc.R") # FIXME remove line when integrated as function in weighted cluster 
		xx <- clueboot(diss, base.clust = base.clust, 
					   k = j, R = R,
					   int.weights = TRUE) 
		
		agg.method <- agg.method
		consAgg <- numeric(length(agg.method))
		for(l in seq_along(agg.method) ){
		  consAgg[l] <- mean(as.vector(cl_agreement(xx, method = agg.method[l])))
		}
		#consAgg <- lapply(consAgg, unlist)
		names(consAgg) <- paste0("cons_",agg.method)
		if(k.fixed){ 
		  cons <- cl_consensus(xx, method = cons.method, control = list(k = j)) 
		} else {
		  cons <- cl_consensus(xx, method = cons.method)
		}
		p()
		if(membership == "fuzzy"){
			list(xx = xx, consAgg = consAgg, clust_ids = cons[[1]])           
		} else {
		  list(xx = xx, consAgg = consAgg, clust_ids = cl_class_ids(cons))           
		} 
	}
  
  clusteringsCons <- lapply(ff, function(x) x$xx)
  consAgreement <- lapply(ff, function(x) x$consAgg)
  consAgreement <- as.data.frame(apply(do.call(rbind, consAgreement), 2, unlist))
  
  
  clustering <- lapply(ff, function(x) x$clust_ids)
  
  names(clustering) <- paste0("cluster", kvals)
  if(membership != "fuzzy"){
  clustering <-  as.data.frame(do.call(cbind, clustering))
  rownames(clustering) <- seq(1, nrow(clustering),1)
  }
  clustResult <- list()
  clustResult$clustering <- clustering
  
  if(membership != "fuzzy"){
    stats <- list()
    
    statsError <-c(0,0,0,0,0,0,0,0,0,1) 
    names(statsError) <-  c("PBC", "HG", "HGSD", "ASW", "ASWw", 
                            "CH", "R2", "CHsq", "R2sq", "HC")
    
    for( i in 1 : ncol(clustering)){
      if(length(unique(clustering[,i])) < 2){
        stats[[i]] <- statsError
        
      }
      if(length(unique(clustering[,i])) >= 2){
        stats[[i]] <- WeightedCluster::wcClusterQuality(diss = diss, 
                                                        clustering = clustering[,i])$stats
      }
    } 
    
    
    stats <- as.data.frame(do.call(rbind, stats), 
                           row.names =  paste0("cluster", kvals) )
    
    stats <- cbind(stats, consAgreement)
    clustResult$stats <- stats  
  }
  
  clustResult$kvals <- kvals
  tStop <- Sys.time()
  dur <- tStop - tStart
  clustResult$call <- sys.call()
  if(keep.ensemble){ 
    clustResult$ensemblePartitions <- clusteringsCons
    names(clustResult$ensemblePartitions) <- paste0("cluster", kvals) 
  }
  cat(paste0("[>] Elapsed time: ", round(dur, digits = 2), 
             " ",attr(dur,"units")), sep = "\n")
  class(clustResult) <- c("consClust")
  return(clustResult)
}

plot.consClust <- function(x, col = NULL, ...){
  if (is.null(col)) {
    allnames <- colnames(x$stats)
    cols <- viridis::turbo(length(allnames))
    names(cols) <- allnames
    cols["RHC"] <- cols["HC"]
    cols <- cols[colnames(x$stats)]
  }
  else{
    cols <- col
  }
  plot.clustrange(x, col = cols, ...)
  
}

print.consClust <- function (x, digits = 2, bootstat = c("t0", "mean", "stderr"), 
          ...) 
{
    x <- round(x$stats, digits)

  print(x, ...)
}

