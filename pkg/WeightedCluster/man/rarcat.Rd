\name{rarcat}
\alias{rarcat}
\title{
Robustness Assessment of Regressions using Cluster Analysis Typologies (RARCAT)
}
\description{
\code{rarcat} is a wrapper for the functions \code{regressboot} and \code{bootpool} that performs the entire RARCAT procedure on all possible associations between a typology and covariates of interest. See Roth et al. (2024) or the R tutorial as \code{WeightedCluster} vignette for all details on the corresponding methods and their utility.
}
\usage{
rarcat(formula, data, diss, 
        robust=TRUE, B=500, count=FALSE, 
        algo="pam", method="ward.D", 
        fixed=FALSE, kcluster=10, cqi="CH",
        parallel="no", ncpus=1, cl=NULL,
        fisher_transform=FALSE, conflevel=0.05, digits=3)
}
\arguments{
  \item{formula}{A formula object with the clustering solution on the left side and the covariates of interest on the ride side.}
  \item{data}{The dataset (data frame) with column names corresponding to the information in formula. The number of individuals (row number) should match the dimension of \code{diss}.}
  \item{diss}{The numerical dissimilarity matrix used for clustering. Only a pre-computed matrix (i.e., where pairwise dissimilarities do not depend on the resample) is currently supported.}
  \item{robust}{Logical. TRUE (the default) indicates that RARCAT should be performed. FALSE implies a much faster function run but only output the original analysis, which is a standard regression analysis for all combinations of reference clusters and covariates.}
  \item{B}{The integer number of bootstrap. Set to 500 by default to attain a satisfactory precision around the estimates as the procedure involves multiple steps.}
  \item{count}{Logical. Whether the bootstrap runs are counted on the screen or not.}
  \item{algo}{The clustering algorithm as a character string. Currently only "pam" (calling the function \code{wcKMedRange}) and "hierarchical" (calling the function \code{fastcluster::hclust}) are supported. By default "pam".}
  \item{method}{A character string with the method argument of \code{hclust}, "ward.D" by default.}
  \item{fixed}{Logical. TRUE implies that the number of clusters is the same in every bootstrap. FALSE (default) implies that an optimal number of clusters is evaluated each time.}
  \item{kcluster}{Integer. Either the number of clusters in every bootstrap if \code{fixed} is TRUE or the maximum number of clusters (starting from 2) to be evaluated in each bootstrap if \code{fixed} is FALSE.}
  \item{cqi}{A character string with the cluster quality index to be evaluated for each new partition. Any column of \code{as.clustrange} is supported, "CH" (the Calinski-Harabasz index) by default. Also works with \code{algo}= "pam".}
  \item{parallel}{A character string with the type of parallel operation to be used (if any) by the function \code{boot:boot}. Options are "no" (default), "multicore" and "snow" (for Windows).}
  \item{ncpus}{Integer. Number of processes to be used in case of parallel operation. Typically, one would chose this to be the number of available CPUs.}
  \item{cl}{A parallel cluster for use if \code{parallel} = "snow". If not supplied, a cluster on the local machine is created for the duration of the \code{boot} call.}
  \item{fisher_transform}{Logical. TRUE means that a Fisher transformation is applied in the \code{bootpool} function. This can be recommended in case of extreme associations (close to the -1 or 1 boundaries). FALSE by default.}
  \item{conflevel}{Confidence level for the confidence intervals from the original analysis and the prediction intervals from the robustness assessment. 0.05 by default.}
  \item{digits}{Controls the number of significant digits to print. 3 by default.}
}
\details{
The \code{rarcat} function runs a standard typology-based association study and evaluates the impact of sampling uncertainty on the results, thus assessing the reproducibility of the analysis.
}
\value{
The output of \code{rarcattables} contains the following tables:
  \item{original.analysis}{Average Marginal Effects (AMEs) estimated with multivariable logistic regressions and representing the expected change in the probability of belonging to a trajectory group (a reference cluster) for a change in the level of a variable (a covariate of interest), together with their confidence intervals.}
  \item{robust.analysis}{Pooled AMEs from the bootstrap procedure and their prediction intervals, representing the range of expected values if the clustering and associated regressions were performed on a new sample from the same underlying distribution. This table provide robust estimates for a typology-based association study.}
}
\references{
Roth, L., Studer, M., Zuercher, E., & Peytremann-Bridevaux, I. (2024). Robustness assessment of regressions using cluster analysis typologies: a bootstrap procedure with application in state sequence analysis. BMC medical research methodology, 24(1), 303. https://doi.org/10.1186/s12874-024-02435-8.
}
\author{
Leonard Roth
}
\seealso{
\code{\link{regressboot}}, \code{\link{bootpool}}
}
\examples{

## Set the seed for reproducible results
set.seed(1)

## Loading the data (TraMineR package)
data(mvad)

## Creating the state sequence object
mvad.seq <- seqdef(mvad, 17:86)

## Distance computation
diss <- seqdist(mvad.seq, method="LCS")

## Hierarchical clustering
hc <- fastcluster::hclust(as.dist(diss), method="ward.D")

## Computing cluster quality measures
clustqual <- as.clustrange(hc, diss=diss, ncluster=6)

## A six clusters solution is chosen here
mvad$clustering <- clustqual$clustering$cluster6

## A formula object with the the covariates of interest (to be related to the typology)
formula <- clustering ~ funemp + gcse5eq

## As in the original analysis, hierarchical clustering with Ward method is implemented
## The number of clusters is fixed to 6 here
## For illustration purposes, the number of bootstrap is smaller than what it ought to be
rarcatout <- rarcat(formula, mvad, diss, B = 50, 
                    algo = "hierarchical", method = "ward.D", 
                    fixed = TRUE, kcluster = 6)

## Assess the robustness of the original analysis
rarcatout$original.analysis
rarcatout$robust.analysis
}
