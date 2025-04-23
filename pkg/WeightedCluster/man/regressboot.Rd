\name{regressboot}
\alias{regressboot}
\title{
Bootstrap replicates of a typology and its association with covariates of interest
}
\description{
The \code{regressboot} function corresponds to the first part of the Robustness Assessment of Regressions using Cluster Analysis Typologies (RARCAT) procedure, which allows for evaluating the impact of sampling uncertainty on a standard Sequence Analysis, and thus assessing the reliability of its findings. See Roth et al. (2024) or the R tutorial as \code{WeightedCluster} vignette for all details on this procedure and its utility. \code{regressboot} should be used together with the \code{unirarcat} function.
}
\usage{
regressboot(diss, covar, df, B = 500, count = FALSE,
            algo = "pam", method = "ward.D",  
            fixed = FALSE, ncluster = 10, eval="CH",
            parallel = "no", ncpus = 1, cl = NULL)
}
\arguments{
  \item{diss}{The numerical dissimilarity matrix used for clustering. Only a pre-computed matrix (i.e., where pairwise dissimilarities do not depend on the resample) is currently supported.}
  \item{covar}{A character vector containing the names of the covariates whose association with the clustering is studied. A formula object is then created inside the function based on this.}
  \item{df}{The dataset (data frame) with the covariates of interest. Column names should match the information in \code{covar} and the number of individuals (row number) should match the dimension of \code{diss}.}
  \item{B}{The integer number of bootstrap. Set to 500 by default to attain a satisfactory precision around the estimates as the procedure involves multiple steps.}
  \item{count}{Logical. Whether the bootstrap runs are counted on the screen or not.}
  \item{algo}{The clustering algorithm as a character string. Currently only "pam" (calling the function \code{wcKMedRange}) and "hierarchical" (calling the function \code{fastcluster::hclust}) are supported. By default "pam".}
  \item{method}{A character string with the method argument of \code{hclust}, "ward.D" by default.}
  \item{fixed}{Logical. TRUE implies that the number of clusters is the same in every bootstrap. FALSE (default) implies that an optimal number of clusters is evaluated each time.}
  \item{ncluster}{Integer. Either the number of clusters in every bootstrap if \code{fixed} is TRUE or the maximum number of clusters (starting from 2) to be evaluated in each bootstrap if \code{fixed} is FALSE.}
  \item{eval}{A character string with the cluster quality index to be evaluated for each new partition. Any column of \code{as.clustrange} is supported, "CH" (the Calinski-Harabasz index) by default. Also works with \code{algo}= "pam".}
  \item{parallel}{A character string with the type of parallel operation to be used (if any) by the function \code{boot:boot}. Options are "no" (default), "multicore" and "snow" (for Windows).}
  \item{ncpus}{Integer. Number of processes to be used in case of parallel operation. Typically, one would chose this to be the number of available CPUs.}
  \item{cl}{A parallel cluster for use if \code{parallel} = "snow". If not supplied, a cluster on the local machine is created for the duration of the \code{boot} call.}
}
\details{
The \code{regressboot} function implements the following steps: (1) A random sample with replacement (i.e, bootstrap) is drawn from the data. (2) The bootstrap sample is clustered applying the exact same clustering procedure as the one used in the original analysis, which implies using the same dissimilarity measure, cluster algorithm, and method to determine the number of clusters. (3) A separate logistic regression predicting membership probability in each group is estimated. (4) The Average Marginal Effect (AME) of each covariate on the probability to be assigned to a given type is retrieved for all sequences belonging to this type. (5) These steps are repeated \code{B} times, with \code{B} typically large.
}
\value{
The output of \code{regressboot} is a list with the following components:
  \item{B}{The number of bootstrap (input parameter).}
  \item{optimal.number}{An integer vector with the numbers of clusters for each bootstrap partition. If input parameter \code{fixed} is FALSE, this corresponds to the selected clustering solution based on the evaluation criterion. If input parameter \code{fixed} is TRUE, this can in rare cases differ from \code{ncluster} if two reference clusters have exactly the same estimated association with a covariate.}
  \item{cluster.solution}{A numerical matrix with the number of individuals (\code{nrow(df)}) as row number and the number of bootstrap (\code{B}) as column number. Each column correspond to the typology for this bootstrap.}
  \item{assoc.char}{A character vector with the different associations evaluated in the logistic regression model (based on input parameter \code{covar}). This corresponds to the name of the covariate for numerical variables and the name with a specific level for factors.}
  \item{original.cluster}{A vector of the same size as the dataset with the original clustering, i.e., the one constructed on the original sample with the given method.}
  \item{original.assoc}{A list with the estimated AMEs corresponding to each association between covariates of interest (as in \code{assoc.char}) and the original typology, i.e., the one constructed on the original sample.}
  \item{coefficients}{A list with the estimated AMEs for all individuals and all bootstraps, corresponding to the associations between covariates of interest (as in \code{assoc.char}) and the typology constructed on each bootstrap. For each covariate, the list contains a numerical matrix with the number of individuals (\code{nrow(df)}) as row number and the number of bootstrap (\code{B}) as column number.}
  \item{errors}{A list with the estimated standard errors of the AMEs for all individuals and all bootstraps, corresponding to the associations between covariates of interest (as in \code{assoc.char}) and the typology constructed on each bootstrap. For each covariate, the list contains a numerical matrix with the number of individuals (\code{nrow(df)}) as row number and the number of bootstrap (\code{B}) as column number.}
}
\references{
Roth, L., Studer, M., Zuercher, E., & Peytremann-Bridevaux, I. (2024). Robustness assessment of regressions using cluster analysis typologies: a bootstrap procedure with application in state sequence analysis. BMC medical research methodology, 24(1), 303. https://doi.org/10.1186/s12874-024-02435-8.

Studer, M. (2013). WeightedCluster library manual: A practical guide to creating typologies of trajectories in the social sciences with R. University of Geneva.

Hennig, C. (2007) Cluster-wise assessment of cluster stability. Computational Statistics and Data Analysis, 52, 258-271.
}
\author{
Leonard Roth
}
\note{
Uses the following packages: fastcluster, dplyr, margins, boot
}
\seealso{
\code{\link{unirarcat}}, \code{\link{rarcat}}
}
\examples{

## Set the seed for reproducible results
set.seed(1)

## Load the margins library for marginal effect estimation
library(margins)

## Loading the data (TraMineR package)
data(mvad)

## Creating the state sequence object
mvad.seq <- seqdef(mvad, 17:86)

## Distance computation
diss <- seqdist(mvad.seq, method="LCS")

## Hierarchical clustering
hc <- fastcluster::hclust(as.dist(diss), method="ward.D")

## Computing cluster quality measures
clustqual <- as.clustrange(hc, diss=diss, ncluster=10)
clustqual

# Create cluster membership variable based on cluster quality above
mvad$clustering <- clustqual$clustering$cluster2
mvad$membership <- mvad$clustering == 2

# Formula for the association between the clustering and a covariate of interest
formula <- membership ~ funemp

# Run logistic regression model
mod <- glm(formula, mvad, family = "binomial")

# Model results
summary(margins(mod))

# A character vector with the name of the covariate of interest (to be related to the typology)
covar <- c("funemp")

## As in the original analysis, hierarchical clustering with Ward method is implemented
## An optimal clustering solution with n between 2 and 10 is evaluated each time by
## maximizing the CH index
## For illustration purposes, the number of bootstrap is smaller than what it ought to be
bootout <- regressboot(diss, covar, mvad, B = 50, 
                      algo = "hierarchical", method = "ward.D", 
                      ncluster = 10)
table(bootout$optimal.number)
bootout$assoc.char
                        
# Robustness assessment for the association between father unemployment status
# and membership to the higher education trajectory group
result <- unirarcat(bootout,  clustqual$clustering$cluster2, 2, "funempyes")
round(result$pooled.ame, 4)
round(result$standard.error, 4)
round(result$bootstrap.deviation, 4)
}
