## function for comparing sets of data by computing LRT and tail probability
##  from the matrix of distances between two groups
## version 1.0, December 2023

# ldiss  matrix of observed dissimilarities or list of such matrices
# group  dichotomous grouping variable or list of grouping variables
# w     vector of case weights or list of such vectors



dissCompare <- function(ldiss,
  lgroup=NULL, #set=NULL,
  s=100, seed=36963, stat="all", squared="LRTonly",
  weighted=TRUE, lw=1,
  #opt=NULL,
  BFopt=NULL,
  #method,
  ...)
{
  #require("gtools")
  #require("TraMineR")

  gc(FALSE)
  ptime.begin <- proc.time()

  if (!is.logical(weighted)) {
    if (weighted != 'by.group')
      stop("weighted must be logical or 'by.group'")
    weight.by <- weighted
    weighted <- TRUE
  }
  else {
    weight.by <- 'global'
  }

  if (is.logical(squared))
    LRTpow <- 1
  else {
    if (squared != "LRTonly") stop("squared must be logical or 'LRTonly'")
    LRTpow <- 2
    squared <- FALSE
  }

  #is1.stslist <- inherits(seqdata,"stslist")
  #is2.stslist <- inherits(seqdata2,"stslist")

  if (!is.list(ldiss)) {
    ldiss <- list(ldiss)
  }
  if (!is.list(lgroup)) {
    lgroup <- list(lgroup)
  }
  if (!is.list(lw))
    lw <- list(lw)

  if (any(!stat %in% c("LRT","BIC","all")))
    stop("Bad stat value, must be one of 'LRT', 'BIC', or 'all'")

  if (any(stat=="all")) {
    is.LRT <- is.BIC <- TRUE
  }
  else{
    is.LRT <- "LRT" %in% stat
    is.BIC <- "BIC" %in% stat
  }

    gvar <- lgroup
    idx.a <- idx.b <- list()
        #seq2 <- list()
    for (i in 1:length(lgroup))


  # prepare samples
  G = length(ldiss)  ## number of sets
  n = matrix(NA,nrow=G,ncol=2)
  #idx.a = idx.b <- list(rep(NA,G))

  for (i in 1:G) {
    if (length(idx.a[[i]])>=length(idx,b[[i]])) {
      n[i,1] <- length(idx.a[[i]])
      n[i,2] <- length(idx,b[[i]])
    }
    else {
      n[i,1] <- length(idx.b[[i]])
      n[i,2] <- length(idx.a[[i]])
    }
  }

  if (s>0) { # for s=0 we do not need that
    m.n = apply(n,1,max)
    n.n = apply(n,1,min)
    f.n1 <- floor(s/m.n)
    ff.n1 <- sapply(f.n1, g<-function(x){max(1,x)})
    #r.n1 = s-m.n%%s
    r.n1 = ifelse(s<m.n, s - m.n%%s, s - f.n1*m.n)
    #k.n = floor((m.n+r.n1)/n.n)
    #r.n2 = (m.n+r.n1)-k.n*n.n
    k.n = floor((ff.n1*m.n+r.n1)/n.n)
    r.n2 = (ff.n1*m.n+r.n1)-k.n*n.n

    ##GR we have an error when s > min(m.n), because then we get some r.n1 > m.n
    ##GR Has been fixed, so we should never reach the following stops
    if(any(m.n<r.n1)) {
      ii <- which(m.n<r.n1)
      stop("rest r.n1 values greater than max m.n for i= ", ii, " s= ", s)
    }
    if(any(n.n<r.n2)) {
      ii <- which(n.n<r.n2)
      stop("rest r.n2 values greater than min n.n for i= ", ii, " s= ", s)
    }
  }

  nc <- ifelse(is.LRT & is.BIC, 4, 2)
  Results <- matrix(NA,G,nc)
  #oopt <- opt
  multsple <- FALSE

  ## permute function from gtools that has become orphaned
  permute <- function(x) sample(x, size = length(x), replace = FALSE)

  ## Constructing vector of indexes of sampled cases
  #r.s1=r.s2 = list(rep(NA,G))
  for (i in 1:G) {
    diss <- ldiss[[i]][c(idx.a[[i]],idx.b[[i]]),c(idx.a[[i]],idx.b[[i]])]
    weights <- lw[[i]][c(idx.a[[i]],idx.b[[i]])]
    if (s==0) { # no sampling
      r1 <- 1:length(idx.a[[i]])
      r2 <- 1:length(idx.b[[i]]) + length(idx.a[[i]])
      #suppressMessages(diss <- seqdist(rbind(seq.a[[i]],seq.b[[i]]), method=method, weighted=weighted, ...))
      #weights <- c(attr(seq.a[[i]],"weights"),attr(seq.b[[i]],"weights"))
    suppressMessages(
        Results[i,] <-
          seq.comp(r1, r2, diss, weights, is.LRT=is.LRT, is.BIC=is.BIC,
             squared=squared, weighted=weighted, weight.by=weight.by,
             LRTpow=LRTpow, ...))
    }
    else { # sampling
      set.seed(seed)
      r.s1 <- c(permute(rep(1:m.n[i],ff.n1[i])),sample(1:m.n[i],r.n1[i],F))
      r.s2 <- c(permute(rep(1:n.n[i],k.n[i])),sample(1:n.n[i],r.n2[i],F))
      r.s1 = matrix(r.s1,ncol=s)
      r.s2 = matrix(r.s2,ncol=s)

      #if (is.null(oopt))
      #  opt <- ifelse(nrow(idx.a[[i]]) + nrow(idx.b[[i]]) > 2*s, 1, 2)
      ##message('opt = ',opt)
      #if (opt==2) {
      #  suppressMessages(diss <- seqdist(rbind(seq.a[[i]],seq.b[[i]]), method=method, weighted=weighted, ...))
      #  weights <- c(attr(seq.a[[i]],"weights"),attr(seq.b[[i]],"weights"))
      #}

      multsple <- nrow(r.s1) > 1 || multsple
    ### new complete samples without replacement of length s over G comparisons
      t<-matrix(NA,nrow=nrow(r.s1),ncol=nc)
      for (j in 1:nrow(r.s1)) {

        #if (opt==2) {
          r1 <- r.s1[j,]
          r2 <- r.s2[j,] + length(idx.a[[i]])
        #}
        suppressMessages(t[j,] <-
            seq.comp(r1, r2, diss, weights, is.LRT=is.LRT, is.BIC=is.BIC,
              squared=squared, weighted=weighted, weight.by=weight.by,
              LRTpow=LRTpow, ...))
      }
      Results[i,]<-apply(t,2,mean)
    }
  }
  colnames <- NULL
  if (is.LRT) colnames <- c("LRT", "p-value")
  if (is.BIC) {
    if (is.null(BFopt) && multsple) {
      BF2 <- exp(Results[,nc-1]/2)
      Results <- cbind(Results, BF2)
      colnames <- c(colnames, "BIC diff.", "Bayes Factor (Avg)", "Bayes Factor (From Avg BIC)")
    }
    else if (BFopt==1 && multsple) {
      colnames <- c(colnames, "BIC diff.", "Bayes Factor (Avg)")
    }
    else if (BFopt==2 && multsple) {
      BF2 <- exp(Results[,nc-1]/2)
      Results[,nc] <- BF2
      colnames <- c(colnames, "BIC diff.", "Bayes Factor (From Avg BIC)")
    }
    else {
      colnames <- c(colnames, "BIC diff.", "Bayes Factor")
    }
  }
  colnames(Results) <- colnames
  if(!is.null(set)) rownames(Results) <- lev.set

  #### Display elaspsed time ####

  ptime.end <- proc.time()
  time.begin <- as.POSIXct(sum(ptime.begin[1:2]), origin = "1960-01-01")
  time.end <- as.POSIXct(sum(ptime.end[1:2]), origin = "1960-01-01")
  time.elapsed <- format(round(difftime(time.end, time.begin), 3))

  message("elapsed time:", time.elapsed)

  return(Results)
}

####################
seq.comp <- function(r1, r2, diss, weights, is.LRT,is.BIC, squared, weighted, weight.by,
                    LRTpow,...)
{
  #print(length(r1))
  #print(length(r2))
  #print(dim(diss))
  #print(length(weights))

  # compute some basic statistics
  #n1 = nrow(S1)
  #n2 = nrow(S2)
  n1 <- length(r1)
  n2 <- length(r2)
  n0 <- n1+n2
  #n.0=log((n0*n0-n0)/2) #not used
  dist.S=dist.S1=dist.S2<-vector()

  weighted <- weighted && !any(is.null(weights))

  ## weights
  if (weighted) {
    w1 <- weights[r1]
    w2 <- weights[r2]
    w <- c(w1,w2)
    ## normalize weights to respect group sizes
    if (weight.by == 'by.group') {
      w1 <- n1/sum(w1) * w1
      w2 <- n2/sum(w2) * w2
      w <- c(w1,w2)
    }
    nw <- sum(w)
    nw1 <- sum(w1)
    nw2 <- sum(w2)
  }
  else { # no weight
    nw <- n0
    nw1 <- n1
    nw2 <- n2
    w <- rep(1,nw)
    w1 <- rep(1,nw1)
    w2 <- rep(1,nw2)
  }

  ## SS <- nw*dissvar(diss[c(r1,r2),c(r1,r2)], weights=w, squared=squared) #
  ## SS1 <- nw1*dissvar(diss[r1,r1], weights=w1, squared=squared) #
  ## SS2 <- nw2*dissvar(diss[r2,r2], weights=w2, squared=squared) #

  ## Using dissvar does not allow the 'LRTonly' solution

  dist.S <-disscenter(diss[c(r1,r2),c(r1,r2)], weights=w, squared=squared) # calculate S distance to center
  dist.S1<-disscenter(diss[r1,r1], weights=w1, squared=squared) # calculate S1 distance to center
  dist.S2<-disscenter(diss[r2,r2], weights=w2, squared=squared) # calculate S2 distance to center

  SS <- sum(w*dist.S^LRTpow)
  SS1 <- sum(w1*dist.S1^LRTpow)
  SS2 <- sum(w2*dist.S2^LRTpow)

  res <- NULL

  LRT <- n0*(log(SS/n0) - log((SS1+SS2)/n0))

  if (is.LRT) {
    p.LRT <- pchisq(LRT,1,lower.tail=FALSE)
    res <- cbind(LRT, p.LRT)
  }

  if (is.BIC) {

    # compute BICs and adjusted BICs
    BIC <- LRT - 1*log(n0)
    # compute Bayes factors
    BF <- exp(BIC/2)

    res <- cbind(res, BIC, BF)
  }
  return(res)
}
