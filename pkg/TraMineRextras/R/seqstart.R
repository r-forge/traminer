seqstart <- function(seqdata, data.start, new.start, tmin=NULL, tmax=NULL, missing=NA){
	new.index <- as.integer(data.start - new.start+1)
	if(length(new.index)!=nrow(seqdata)){
		if(length(new.index)==1){
			new.index <- rep(new.index, nrow(seqdata))
		}
		else{
			stop(" [!] number of individual in seqdata, data.start and/or new.start mismatch.")
		}
	}
	if(any(new.index <1)){
		correction <- -min(new.index)+1
		new.index <- new.index + correction
	} else {
		correction <- 0
	}
	if(is.null(tmin)){
		tmin <- max(min(new.index), 1)
	}
	tmin <- tmin+correction
	if(is.null(tmax)){
		tmax <- max(new.index)+ncol(seqdata) -1
	}
	tmax <- tmax +correction
	if(tmax<0){
		stop("[!] There are no data in the specified new time frame.")
	}
	## cat(tmin, tmax, correction)
	##new.index.mat <- new.index - tmin + 1
	new.data <- matrix(as.character(missing), ncol=(tmax-tmin+1), nrow=nrow(seqdata))
	seqindex <- matrix(tmin:tmax, ncol=(tmax-tmin+1), nrow=nrow(seqdata), byrow=TRUE)
	seqdataindex <- matrix(1:ncol(seqdata), ncol=ncol(seqdata), nrow=nrow(seqdata), byrow=TRUE)
	seqdata <- as.matrix(seqdata)
	
	rowindex <- (1:ncol(seqdata))-1
	for(ind in 1:nrow(seqdata)){
		indexes <-  new.index[ind]+rowindex - tmin +1
		cond <- (indexes <= ncol(new.data)) & (indexes > 0)
		## print(indexes)
		new.data[ind, indexes[cond]] <- seqdata[ind, cond]
		#print(seqdata[which.indiv, oldindex])
		##print(new.data)
	}
	return(new.data)
}

