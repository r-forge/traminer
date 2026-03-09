	\name{consClust}
\alias{consClust}


\title{
	Crisp Consensus Clustering 
}
\description{

Compute consensus clustering for different number of clusters (Monti, 2003). The function further computes cluster quality and consensus agreement measures. 

}

\usage{
consClust(diss,
          base.clust = "pam", 
          R = 100, 
          kvals = 2:15,
          cons.method = "SE", 
          membership = "crisp",
          k.fixed = TRUE, 
          agg.method = "cRand",
          keep.ensemble = TRUE,
          parallel = FALSE,
          progressbar = TRUE)
}

\arguments{
  \item{diss}{A dissimilarity matrix or a \code{dist} object.}
  
  \item{base.clust}{Character. Clustering algorithms used to compute the ensemble of partitions and hierarchies. May be a combination of \code{"pam", "single", "complete", "average", "mcquitty", "ward.D", "ward.D2", "centroid", "median"}.}
  
  \item{R}{Numeric. The number of partitions or hierarchies to compute for consensus clustering.}
  
  \item{kvals}{Numeric vector. The number of clusters to compute, default \code{2:15}}
  
  \item{cons.method}{Character. The consensus clustering method to use, can be one of \code{"SE"} (default), \code{"HE", "SM", "HM", "GV1", "DWH", "GV3", "soft/symdiff", "hard/symdiff"}. See \code{\link[clue]{cl_consensus}} for details on the methods.}
  
  \item{membership}{Character. If \code{"crisp"}, the consensus clustering is returned as vectors of crisp cluster labels. If \code{"fuzzy"} the function returns fuzzy membership martices.}
  
  \item{k.fixed}{Logical. If \code{TRUE} (default), the number of clusters obtained from the consensus cannot exceed the number of cluster in the partition ensemble.}
  
  \item{agg.method}{Character. The consensus agreement measures to compute, may be a comination of \code{"cRand"} (default), \code{"Rand", "euclidean", "manhattan", "NMI", "KP", "angle", "diag", "FM", "Jaccard", "puritiy", "PS"}. See \code{\link[clue]{cl_agreement}} for details on the methods.}
  
  \item{keep.ensemble}{Logical. If \code{TRUE} (default) partitions and/or hierarchies are returned by the function. Setting \code{keep.ensemble = FALSE} saves memory.}
  
  \item{parallel}{Logical. Whether to initialize the parallel processing of the \code{future} package using the default \code{\link[future]{multisession}} strategy. If \code{FALSE} (default), then the current \code{\link[future]{plan}} is used. If \code{TRUE}, \code{\link[future]{multisession}} \code{\link[future]{plan}} is initialized using default values.}
  
	\item{progressbar}{Logical. Whether to initialize a progress bar using the \code{future} package. If \code{FALSE} (default), then the current progress bar \code{\link[progressr]{handlers}} is used . If \code{TRUE}, a new global progress bar \code{\link[progressr]{handlers}} is initialized.}
	}

\details{
\code{consClust} relies on \code{\link[clue]{cl_consensus}}, to compute a consensus clustering among several internally computed hierarchies and partitions. The algorithm works as follows: 

\enumerate{
  \item An ensemble of \code{R} clusterings in a fixed number of groups \code{k} in \code{kvals} are computed on subsamples of the data. To reflect the potential data perculiarities, clustering are obtained on weighted subsamples using Baysian resampling.
  
  \item A consensus among the clusterings obtained in step one is computed. The number of clusters in the consensus may exceed the one in the ensemble of clusterings computed in step one. Setting \code{k.fixed} to \code{TRUE} set the maximal number of cluster in the consensus to the number of clusters \code{k} in the ensemble.
  
  \item Cluster quality indices are computed for the obtained consensus.
  
  \item Step 1 to 3 are repeated for each number of groups specified in \code{kvals}
}


}
\value{
	A \code{consClust} object with the following components:
	\item{\code{clustering}}{The retained clustering for each number of groups.}
	\item{\code{stats}}{A \code{matrix} containing the clustering statistics of each cluster solution.}
	\item{\code{kvals}}{The number of computed clusters.}
	\item{\code{call}}{The used function calls.}
	\item{\code{ensemblePartitions}}{A list containing the partitions or hierarchies used to obtain the consensus.}
}

\references{
Monti, S., Tamayo, P., Mesirov, J., Golub, T. (2003). Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data. \emph{Machine Learning, 52, 1} \doi{10.1023/A:1023949509487}
  
Unterlerchner, L., Studer, M. (2026). What are We Looking For? A Comparative Review of Clustering Algorithms and Cluster Quality Indices For Sequence Analysis. \emph{LIVES Working Papers 108} \doi{10.12682/lives.2296-1658.2026.108}
	
}

\examples{
# Loading illustrative data
data(mvad)

# Creating state sequence object
mvad.seq <- seqdef(mvad[1:200, 17:86])

# Computing dissimilarities using LCS measure
diss <- seqdist(mvad.seq, method="LCS")

## Computing consensus clustering using PAM and Ward (D)

pamWardConsClust <- consClust(diss,
                              kvals = 2:6, 
                              base.clust = c("pam", "ward.D"),
                              R = 20,
                              k.fixed = TRUE,
                              agg.method = "cRand")

## Showing the cluster quality measures. 
pamWardConsClust

## Plotting normalized values for easier identification of minimum and maximum values, with a transparent legend background.
plot(pamWardConsClust, norm="range")


# Plotting sequences in 6 groups
par(mar = c(2.5,2,1.8,1.2))

seqdplot(mvad.seq, 
         group = pamWardConsClust$clustering$cluster6, 
         border = NA)
}
