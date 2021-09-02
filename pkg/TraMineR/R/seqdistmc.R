## multichannel distances

seqdistmc <- function(channels, method, norm="none", indel=1, sm=NULL,
	with.missing=FALSE, full.matrix=TRUE, link="sum", cval=2, miss.cost=2, cweight=NULL,
  what="diss", ch.sep = "@@@@TraMineRSep@@@@") {

	## Checking arguments
  whatlist <- c("diss","sm","seqmc")
  if (!(what %in% whatlist)){
    msg.stop("what should be one of ",paste0("'",whatlist,"'", collapse=","))
  }

  if(length(indel) > 1 & any(indel=="auto"))
    stop(" [!] 'auto' not allowed in vector or list indel")

	nchannels <- length(channels)
	if (nchannels < 2) {
		stop("[!] please specify at least two channels")
	}
	if (is.null(cweight)) {
		cweight <- rep(1, nchannels)
	}
  for (i in 1:nchannels){
    if (length(grep(ch.sep, alphabet(channels[[i]], with.missing=TRUE), fixed=TRUE))>0)
      stop(" [!] ch.sep symbol (",ch.sep,") occurs in alphabet of at least one channel")
  }
  if (is.list(indel) & length(indel) != nchannels)
		stop("[!] when a list, indel must be of length equal to number of channels")
	numseq <- sapply(channels,nrow)
	if(any(numseq!=numseq[1])) {
		stop(" [!] sequence objects have different numbers of rows")
	}
	numseq <- numseq[1]
	message(" [>] ", nchannels, " channels with ", numseq, " sequences")
	## Actually LCP and RLCP are not included
	metlist <- c("OM", "LCS", "DHD", "HAM")
	if (!method %in% metlist) {
		stop(" [!] method must be one of: ", paste(metlist, collapse=" "), call.=FALSE)
	}
	## We handle LCS as OM
	if (method=="LCS") {
		method <- "OM"
		sm <- "CONSTANT"
		indel <- 1
		cval <- 2
		miss.cost <- 2
	}
	timeVarying <- method %in% c("DHD")
	## Indels and sm only apply for OM
	## Correct number of arguments
	## Generating default substitution arguments for DHD and HAM
	if (is.null(sm)) {
		costmethod <- "CONSTANT"
		if (method == "DHD") {
			costmethod <- "TRATE"
		}
		sm <- rep(costmethod, nchannels)
	}
	else if (length(sm)==1 && sm %in% c("CONSTANT", "TRATE", "INDELS", "INDELSLOG")){
		sm <- rep(sm, nchannels)
	}

	if (length(indel)==1) {
	   	indel <- rep(indel, nchannels)
	}
	## Checking correct numbers of info per channel
	if ((length(indel)!= nchannels) ||
		(length(sm)!= nchannels) ||
		(length(cweight)!= nchannels)) {
		stop(" [!] you should supply one weight, one substitution matrix and one indel per channel")
	}
	## indels
	if (!is.list(indel))
    indel_list <- numeric(length=nchannels)
  if (any(indel == "auto") & any(sm %in% c("INDELSLOG","INDELS")))
    indel_list <- list()
	## subsitution matrix
	substmat_list <- list()
	## alphabet for each channel
	alphabet_list <- list()
	## alphabet size per channel
	alphsize_list <-list()
	## seqlenth of each channels
	maxlength_list <- numeric(length=nchannels)

  if (what != "seqmc") {
  	## ============================================================
  	## Building and checking substitution matrix per channel
  	## ============================================================
  	for (i in 1:nchannels) {
  		## Sequence object
  		if (!inherits(channels[[i]],"stslist")) {
  			stop(" [!] channel ", i, " is not a state sequence object, use 'seqdef' function to create one", call.=FALSE)
  		}
  		alphabet_list[[i]] <- attr(channels[[i]],"alphabet")
  		## Checking missing values
  		if (with.missing) {
  			alphabet_list[[i]] <- c(alphabet_list[[i]],attr(channels[[i]],"nr"))
  			message(" [>] including missing value as an additional state" )
  		}
  		else {
  			if (any(channels[[i]]==attr(channels[[i]],"nr"))) {
  				stop(" [!] found missing values in channel ", i, ", please set 'with.missing=T' to nevertheless compute distances")
  			}
  		}
  		alphsize_list[[i]] <- length(alphabet_list[[i]])
      if(is.list(indel)){
        if (length(indel[[i]])==1)
          indel[[i]] <- rep(indel[[i]],alphsize_list[[i]])
        if (length(indel[[i]]) != alphsize_list[[i]])
  				stop(" [!] indel length does not much size of alphabet for at least one channel")
      }
      else if (!any(indel=="auto") & !is.list(indel_list)) {
  		  indel_list[i] <- indel[i]
      }

  		## Storing number of columns
  		maxlength_list[i] <- ncol(channels[[i]])
  		## Substitution matrix generation method is given
  		if	(is.character(sm[[i]])) {
  			message(" [>] computing substitution cost matrix for channel ", i)
  			costs <- seqcost(channels[[i]], sm[[i]], with.missing=with.missing,
  				time.varying=timeVarying, cval=cval, miss.cost=miss.cost)
        substmat_list[[i]] <- costs$sm
        if (any(indel=="auto")) {
          if (is.list(indel_list))
            indel_list[[i]] <- costs$indel
          else
            indel_list[i] <- costs$indel
        }
  		}
  		## Checking correct dimension cost matrix
  		else {
  			if (method=="OM") {
          if (any(indel[i] == "auto"))
              indel_list[i] <- max(sm[[i]])/2
          else
              indel_list[i] <- indel[i]
          cat("\n indel_list[i] ",indel_list[i], "\n")
          print(sm[[i]])
  				checkcost(sm[[i]], channels[[i]], with.missing = with.missing, indel = indel_list[i])
  			} else {
  				checkcost(sm[[i]], channels[[i]], with.missing = with.missing)
  			}
  			substmat_list[[i]] <- sm[[i]]
  		}

  		## Mutliply by channel weight
  		substmat_list[[i]] <- cweight[i]* substmat_list[[i]]
  	}
    if (any(indel=="auto")) indel <- indel_list
  } else {
      for (i in 1:nchannels) {
  		  maxlength_list[i] <- ncol(channels[[i]])
      }
  }

	## Checking that all channels have the same length
	slength1 <- seqlength(channels[[1]])
	for (i in 2:nchannels) {
		if (sum(slength1 != seqlength(channels[[i]]))>0) {
			if (!with.missing) {
				stop(" [!] Some individuals have channels of different length. Set 'with.missing=TRUE'.")
			} else {
				warning(" [!] Some individuals have channels of different length. Shorter sequences will be filled with missing values.")
				break
			}
		}
	}
	## ================================
	## Building the new sequence object
	## ================================
	message(" [>] building combined sequences...", appendLF=F)
	## Complex separator to ensure (hahem) unicity
	##sep <- "@@@@TraMineRSep@@@@"  ## now argument ch.sep
  sep <- ch.sep
	maxlength=max(maxlength_list)
	newseqdata <- matrix("", nrow=numseq, ncol=maxlength)
  rownames(newseqdata) <- rownames(channels[[1]])
	newseqdataNA <- matrix(TRUE, nrow=numseq, ncol=maxlength)
	for (i in 1:nchannels) {
		seqchan <- channels[[i]]
		void <- attr(seqchan, "void")
		nr <- attr(seqchan, "nr")
		for (j in 1:maxlength) {
			## No column in stslist object, filling with voids
			if (j > maxlength_list[i]) {
				newCol <- as.character(rep(void, numseq))
			}
			else {
				newCol <- as.character(seqchan[,j])
			}
			## If all channel are equal to void, then we accept as void
			newseqdataNA[,j] <- newseqdataNA[,j] & newCol == void
			## Setting void as nr
			newCol[newCol == void] <- nr
			if (i > 1) {
				newseqdata[,j] <- paste(newseqdata[,j], newCol, sep = sep)
			}
			else {
				newseqdata[,j] <- newCol
			}
		}
  }
	## Setting void states back to NA  (nr will be considered as a distinct state)
	newseqdata[newseqdataNA] <- NA

  ## since v 2.2-0 automatic cpal no longer limited to 12 states, so no need of following
	#alphabet_size <- length(unique(as.character(newseqdata))) - as.integer(sum(is.na(newseqdata))>0)
	#suppressMessages(newseqdata <- seqdef(newseqdata, cpal=rep("blue", alphabet_size)))
	suppressMessages(newseqdata <- seqdef(newseqdata))
	message(" OK")

  if (what == "seqmc") {
    return(newseqdata)
  } else {
  	## =========================================
  	## Building the new substitution cost matrix
  	## =========================================
  	message(" [>] computing combined substitution and indel costs...", appendLF=FALSE)
  	## Build subsitution matrix and new alphabet
  	alphabet <- attr(newseqdata,"alphabet")
  	alphabet_size <- length(alphabet)
    newindel <- NULL
  	## Recomputing the subsitution matrix
  	if (!timeVarying) {
  		newsm <- matrix(0, nrow=alphabet_size, ncol=alphabet_size)
      if (is.list(indel)){
        newindel <- rep(0,alphabet_size)
        statelisti <- strsplit(alphabet[alphabet_size], sep, fixed=TRUE)[[1]]
        for (chan in 1:nchannels){
					 ipos <- match(statelisti[chan], alphabet_list[[chan]])
           newindel[alphabet_size] <- newindel[alphabet_size] + indel[[chan]][ipos]*cweight[chan]
        }
      }
  		for (i in 1:(alphabet_size-1)) {
  			statelisti <- strsplit(alphabet[i], sep, fixed=TRUE)[[1]]
        if (is.list(indel)){
          for (chan in 1:nchannels){
  					 ipos <- match(statelisti[chan], alphabet_list[[chan]])
             newindel[i] <- newindel[i] + indel[[chan]][ipos]*cweight[chan]
          }
        }
  			for (j in (i+1):alphabet_size) {
  				cost <- 0
  				statelistj <- strsplit(alphabet[j], sep, fixed=TRUE)[[1]]
  				for (chan in 1:nchannels) {
  					ipos <- match(statelisti[chan], alphabet_list[[chan]])
  					jpos <- match(statelistj[chan], alphabet_list[[chan]])
  					cost <- cost + substmat_list[[chan]][ipos, jpos]
  				}
  				newsm[i, j] <- cost
  				newsm[j, i] <- cost
  			}
  		}
  	} else {
  		## Recomputing time varying substitution
  		newsm <- array(0, dim=c(alphabet_size, alphabet_size, maxlength))
  		for (t in 1:maxlength) {
  			for (i in 1:(alphabet_size-1)) {
  				statelisti <- strsplit(alphabet[i], sep, fixed=TRUE)[[1]]
  				for (j in (i+1):alphabet_size) {
  					cost <- 0
  					statelistj <- strsplit(alphabet[j], sep, fixed=TRUE)[[1]]
  					for (chan in 1:nchannels) {
  						ipos <- match(statelisti[chan], alphabet_list[[chan]])
  						jpos <- match(statelistj[chan], alphabet_list[[chan]])
  						cost <- cost + substmat_list[[chan]][ipos, jpos, t]
  					}
  					newsm[i, j, t] <- cost
  					newsm[j, i, t] <- cost
  				}
  			}
  		}
  	}
    rownames(newsm) <- colnames(newsm) <- alphabet  ## labels are too long
  	message(" OK")
  	## Indel as sum
    if (is.null(newindel) & !is.list(indel_list)) {
  	   newindel <- sum(indel*cweight)
    }
  	## If we want the mean of cost..
  	if (link=="mean") {
  		newindel <- newindel / sum(cweight)
  		newsm <- newsm / sum(cweight)
  	}
    if (what == "sm") {
      attr(newsm,"indel") <- newindel
      attr(newsm,"alphabet") <- alphabet
      attr(newsm,"cweight") <- cweight
      return(newsm)
    }
  }
  if (what == "diss") {
	   message(" [>] computing distances ...")
	   ## Calling seqdist
	   return(seqdist(newseqdata, method=method, norm=norm, indel=newindel,
		        sm=newsm, with.missing=FALSE, full.matrix=full.matrix))
  }
}
