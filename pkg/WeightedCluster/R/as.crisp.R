as.crisp <- function(clara){
 fuzzClust <- clara
  crispClust <- list()
  labels <- list()
  for(i in 1:length(fuzzClust$clustering)){
    crispClust[[i]] <- apply(fuzzClust$clustering[[i]], 1, which.max)
    
    if(length(unique(crispClust[[i]])) == (fuzzClust$kvals[[i]]+1)){
      labels <- colnames(fuzzClust$clustering[[i]])
    }
      else{
        labels <- colnames(fuzzClust$clustering[[i]])[1:fuzzClust$kvals[[i]]]
      }
     crispClust[[i]] <- as.character(factor(crispClust[[i]], # trouver solution pour ce bug, il y a des cas ou aucune observation est classee dans noise
                                     labels = labels))
  }
  crispClust <- as.data.frame(do.call(cbind, crispClust))
  colnames(crispClust) <- names(fuzzClust$clustering)
  fuzzClust$clustering <- crispClust
  return(fuzzClust)
}

