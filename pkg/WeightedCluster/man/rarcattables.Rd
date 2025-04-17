\name{rarcattables}
\alias{rarcattables}
\title{
Robustness Assessment of Regressions using Cluster Analysis Typologies (RARCAT)
}
\description{
\code{rarcattables} is a wrapper for the function \code{rarcat} that generalizes it to all possible associations between a typology and covariates of interest. See Roth et al. (2024) or the two R tutorials as \code{WeightedCluster} package vignettes for all details on the corresponding methods and their utility. \code{rarcattables} should be used together with the \code{regressboot} function.
}
\usage{
rarcattables(clustering, covar, df, bootout=NULL, 
            transformation=FALSE, conflevel=0.05, digits=3)
}
\arguments{
  \item{clustering}{An integer vector containing the clustering solution (one entry for each individual) from the original analysis.}
  \item{covar}{A character vector containing the names of the covariates whose association with the clustering is studied. The formula object is then created inside the function based on this.}
  \item{df}{The original dataset (data frame) containing the covariates of interest. Row number should be equal to the length of the \code{clustering} argument and column names should match the information in \code{covar}.}
  \item{bootout}{Output of the \code{regressboot} function. Optional, if not given the function simply runs a standard regression analysis for all combinations of reference clusters and covariates.}
  \item{transformation}{Logical. TRUE means that a Fisher transformation is applied in the \code{rarcat} function. This can be recommended in case of extreme associations (close to the -1 or 1 boundaries). FALSE by default.}
  \item{conflevel}{Confidence level for the confidence intervals from the original analysis and the prediction intervals from the robustness assessment. 0.05 by default.}
  \item{digits}{Controls the number of significant digits to print. 3 by default.}
}
\details{
The \code{rarcattables} function runs a standard typology-based association study and evaluates the impact of sampling uncertainty on the results, thus assessing the reproducibility of the analysis.
}
\value{
The output of \code{rarcattables} containing the following tables:
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
\code{\link{regressboot}}, \code{\link{rarcat}}
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

# A character vector with the names of the covariates of interest (to be related to the typology)
covar <- c("funemp", "gcse5eq")

## As in the original analysis, hierarchical clustering with Ward method is implemented
## The number of clusters is fixed to 6 here
## For illustration purposes, the number of bootstrap is smaller than what it ought to be
bootout <- regressboot(diss, covar, mvad, B = 50, 
                        algo = "hierarchical", method = "ward.D", 
                        fixed = TRUE, ncluster = 6)
                        
## Assess the robustness of the original analysis
rarcatout <- rarcattables(clustqual$clustering$cluster6, covar, mvad, bootout)
rarcatout$original.analysis
rarcatout$robust.analysis
}
