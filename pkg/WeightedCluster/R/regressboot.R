
regressboot <- function(formula, data, diss, clustering, factorName, R=500, kmedoid = FALSE, hclust.method="ward.D",
                        fixed=FALSE, ncluster=10, cqi="HC", parallel=FALSE, progressbar=FALSE, debug=FALSE) {
  
  # Ensure dissimilarity matrix and dataset have the same size
  if(nrow(data) != nrow(diss)) {
    stop(" [!] Please verify that the dissimilarity matrix corresponds to the dataset")
  }
  if(length(ncluster) > 1) {
    stop(" [!] Please give a single number as kcluster (see documentation)")
  }
  
  if (parallel) {
    message(" [>] Initializing parallel processing...", appendLF = FALSE)
    oplan <- plan(multisession)
    on.exit(plan(oplan), add = TRUE)
    message(" OK.")
  }
  if (progressbar) {
    if (requireNamespace("progress", quietly = TRUE)) {
      old_handlers <- handlers(handler_progress(format = "(:spin) [:bar] :percent | Elapsed: :elapsed | ETA: :eta | :message"))
      if (!is.null(old_handlers)) {
        on.exit(handlers(old_handlers), add = TRUE)
      }
    } else {
      message(" [>] Install the progress package to see estimated remaining time.")
    }
    oldglobal <- handlers(global = TRUE)
    if (!is.null(oldglobal)) {
      on.exit(handlers(global = oldglobal), add = TRUE)
    }
  }
  p <- progressor(R)

  ## Memory clean up before parallel computing
  gc()
  message(" [>] Starting bootstraps...\n")
  
  ## Random component made at the beginning relies on current seed, no parallel issues.
  ## Should we add the original indices in the first bootstraps ?
  boot_indices <- replicate(R, sample(seq_len(nrow(data)), replace = TRUE), simplify = FALSE)
  #boot_indices[[1]] <- seq_len(nrow(data))
  
	res_list <- foreach(i = seq_len(R), .options.future = list(seed = TRUE, 
						globals=c("diss", "boot_indices", "formula", "data", "kmedoid", "hclust.method",
								"fixed", "ncluster", "cqi", "debug",  "p"),  # reproducible RNG
					packages = c("WeightedCluster", "stats", "progressr"))) %dofuture% {

	  # Subsample indices
	  indices <- boot_indices[[i]]
	  
	  # Bootstrap
	  ame <- amemat(diss = diss,
					indices=indices,
					formula = formula,
					data = data,
					kmedoid = kmedoid, 
					hclust.method=hclust.method,
					fixed = fixed, 
					ncluster=ncluster, 
					cqi=cqi, 
					debug=debug)
		p()
		ame, nrow=nrow(data)
	}
 message(" [>] Finished bootstraps.\n")
  
  #res <- boot::boot(diss, amemat, R, formula=formula, data1=data, 
#					kmedoid = kmedoid, hclust.method=hclust.method, cqi=cqi, ncluster=ncluster, fixed=fixed)
  originres <- res$t0
  bootres <- t(res$t)
  
  cluster.solution <- sapply(res_list, function(x) (x)[, "clustering"]) 
  optimal.number <- apply(cluster.solution, 2, function(x) length(unique(na.omit(x))))
   
  original.assoc <- list()
  effects <- list()
  errors <- list()
  ids <- seq_along(clustering)
  message(" [>] Pooling results...")
  for (ff in factorName) {
    effects[[ff]] <- sapply(res_list, function(x) (x)[, ff]) 
    errors[[ff]] <- sapply(res_list, function(x) (x)[, paste(ff, "SE")])
	for(cc in unique(clustering)){
		ame <- effects[[ff]][clustering == cc,]
		error <- errors[[ff]][clustering == cc,]
		# Prepare the data to input in the multilevel model
		prep <- data.frame(bootstrap = rep(1:ncol(ame), nrow(ame)),
						 id = rep(ids[clustering == cc], each=ncol(ame)),
						 ame = c(t(ame)), 
						 sterror = c(t(error)))
		prep <- subset(prep, !is.na(ame))
		
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
  }
  message("OK.")
  
  
  ## Directly do boot pool
  
  
	  
	  

  }
  ret <- list(list("B" = B, 
              "optimal.kcluster" = optimal.number, 
              "cluster.solution" = cluster.solution, 
              "covar.name" = varlist,
              "original.cluster" = original.cluster,
              "original.ame" = original.assoc,
              "bootstrap.ame" = effects, 
              "std.err" = errors))
	class(ret) <- "regressboot"
	return(ret)
}
