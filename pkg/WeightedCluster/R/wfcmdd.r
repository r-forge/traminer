wfcmdd <- function (diss, memb, weights=NULL, method = "FCMdd", 
    m = 2, dnoise = NULL, eta = NULL, alpha = 0.001, iter.max = 100,  verbose = FALSE) 
{

	## Setting and checking arguments values 
    METHODS <- c("NCdd", "HNCdd", "FCMdd", "PCMdd")
    method <- match.arg(method, METHODS)
    if (method == "KMdd") {
        m = 1
        dnoise = NULL
        eta = NULL
    }
    else if (method == "FCMdd") {
        dnoise = NULL
        eta = NULL
    }
    else if (method == "NCdd") {
        if (is.null(dnoise)) 
            stop("Must provide a value for dnoise")
        eta = NULL
    }
    else if (method == "HNCdd") {
        if (is.null(dnoise)) 
            stop("Must provide a value for dnoise")
        eta = NULL
        m = 1
    }
    else if (method == "PCMdd") {
        if (is.null(eta)) 
            stop("Must provide a vector of values for eta")
        dnoise = NULL
    }
    d = as.matrix(diss)
    n = nrow(d)
    if (is.data.frame(memb) || is.matrix(memb)) {
        if (nrow(memb) != ncol(d)) {
            stop("The number of rows in memb must be the same as the number rows and columns of d")
        }
        u = as.matrix(memb)
    }
    else if (is.vector(memb) && is.numeric(memb)) {
        u = matrix(0, n, length(memb))
        for (k in 1:length(memb)) {
            u[memb[k], k] = 1
        }
    }
    else {
        stop("Provide a number, a vector of seeds, or membership matrix for mobile clusters")
    }
    kMov = ncol(u)
    
    med = rep(NA, kMov)
    if ((method == "PCMdd") && length(eta) != kMov) 
        stop("Vector of reference distances (eta) must have a length equal to the number of clusters")
    if (method == "NCdd" || method == "HNCdd") {
        u = cbind(u, vector("numeric", length = n))
    }
    uPrev = matrix(0, nrow = n, ncol = kMov)
    
		dist2med = matrix(0, nrow = n, ncol = kMov)
    
    continue = TRUE
    iter = 1
    while (continue) {
			for (k in 1:kMov) {
				candidates <- which(apply(u[, -k, drop = FALSE],  1, max) < 1)
			  if(k>1){ ## Avoid previously chosen medoids.
			  	candidates <- candidates[!candidates %in% med[1:(k-1)]]
			  }
			  med[k] = candidates[which.min((u[, k]^m*weights) %*% d[, candidates])]
			  dist2med[, k] = d[, med[k]]
			}
		
			if (method == "KMdd") {
                minC <- apply(dist2med, 1, which.min)
                u[, ] = 0
                for (k in 1:length(minC)) u[k, minC[k]] = 1
            }
            
            else if (method == "HNCdd") {
                d2cm <- cbind(dist2med, dnoise)
                u[, ] = 0
                minC <- apply(d2cm, 1, which.min)
                for (k in 1:length(minC)) {
                  u[k, minC[k]] = 1
                }
            }
            else if (method == "NCdd") {
                d2cm <- cbind(dist2med, dnoise)
                for (k in 1:ncol(d2cm)) {
                  a <- sweep(d2cm, 1, d2cm[, k], "/")
                  u[, k] = 1/rowSums(a^(-1/(m - 1)))
                }
                u[d2cm == 0] = 1
            }
            else if (method == "FCMdd") {
                d2cm <- dist2med
                for (k in 1:ncol(d2cm)) {
                  a <- sweep(d2cm, 1, d2cm[, k], "/")
                  u[, k] = 1/rowSums(a^(-1/(m - 1)))
                }
                u[dist2med == 0] = 1
            }
            else if (method == "PCMdd") {
                for (k in 1:ncol(dist2med)) u[, k] = 1/(1 + ((dist2med[, 
                  k])/eta[k])^(1/(m - 1)))
                u[dist2cent == 0] = 1
            }
            if (iter > 2) {
                continue = (max(abs(u - uPrev)) > alpha) && (iter <= 
                  iter.max) && (max(abs(u - uPrev2)) > alpha)
            }
            if (continue) {
                uPrev2 = uPrev
                uPrev = u
                iter = iter + 1
                if (verbose) 
                  cat(".")
            }
        }
		if (method == "FCMdd") 
            functional = sum(dist2med * (u^m) *weights)
      else if (method == "NCdd" || method == "HNCdd") 
          functional = sum(dist2med * (u[, -(kMov + 
              1)]^m*weights)) + sum(dnoise * u[, kMov +  1]^m *weights)
      else if (method == "PCMdd") {
          functional = 0
          for (k in 1:(kMov + 1)) functional = functional + 
              sum(dist2med[, k] * (u[, k]^m)*weights) + sum(eta[k] * 
              (1 - u[, k])^m*weights)
      }
      if (verbose) 
          cat(paste("\nIterations:", iter, "Functional: ", 
              functional, "\n"))
      u = as.data.frame(u)
      dist2cent = as.data.frame(dist2med)
        
      for (k in 1:kMov) {
          names(u)[k] = paste("M", k, sep = "")
          names(dist2cent)[k] = paste("M", k, sep = "")
      }
   
      if (method == "NCdd" || method == "HNCdd") 
          names(u)[kMov  + 1] = "N"
      rownames(u) = rownames(d)
      rownames(dist2cent) = rownames(d)
      size = colSums(u[, 1:(kMov)])
      withinss = colSums((dist2cent) * (u[, 1:kMov]^m*weights))
      mobileCenters = NULL
      if (method == "KMdd" || method == "FCMdd" || method == 
          "NCdd" || method == "HNCdd" || method == "PCMdd") {
          mobileCenters = med[1:kMov]

      }
      res = list(mode = "dist", method = method, m = m, dnoise = dnoise, 
          eta = eta, memb = u, mobileCenters = mobileCenters, 
          fixedCenters = NULL, dist2clusters = dist2cent, 
          withinss = withinss, size = size, functional = functional)
      class(res) <- "vegclust"
      return(res)

}
