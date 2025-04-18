\name{rarcat}
\alias{rarcat}
\title{
Bootstrap-based robustness assessment for an association between a cluster and a covariate
}
\description{
The \code{rarcat} function corresponds to the second part of the Robustness Assessment of Regressions using Cluster Analysis Typologies (RARCAT) procedure, which allows for evaluating the impact of sampling uncertainty on a standard Sequence Analysis, and thus assessing the reliability of its findings. See Roth et al. (2024) or the two R tutorials as \code{WeightedCluster} package vignettes for all details on this procedure and its utility. \code{rarcat} should be used together with the \code{regressboot} function.}
\usage{
rarcat(bootout, clustering, clusnb, assoc, transformation = FALSE)
}
\arguments{
  \item{bootout}{Output of the \code{regressboot} function.}
  \item{clustering}{An integer vector containing the clustering solution (one entry for each individual) from the original analysis.}
  \item{clusnb}{An integer with the cluster to be evaluated (part of the clustering solution), as the RARCAT procedure is cluster-wise by design.}
  \item{assoc}{A character string with the association of interest as specified in the component \code{assoc.char} of the \code{regressboot} function output.}
  \item{transformation}{Logical. TRUE means that the Average Marginal Effects (AMEs) from the bootstrap procedure are transformed with a Fisher transformation before being imputed in the pooling model, and then transformed back for the output results. This can be recommended in case of extreme associations (close to the -1 or 1 boundaries). FALSE by default.}
}
\details{
The \code{rarcat} function takes as input the AMEs (for each individual and each bootstrap) and their standard errors estimated with the \code{regressboot}. It then combine them using a multilevel modelling framework that mimics a meta-analysis. The summary estimates of effect thus produced account for the sampling uncertainty and should be compared with the results from the original analysis to assess their robustness. Moreover, the individual random effects inform on the central and outlier trajectories in a cluster.
}
\value{
The output of \code{rarcat} is a list with the following components:
  \item{nobs}{An integer with the number of observations (i.e., number of estimated AMES from the function \code{regressboot}) used to compute the robust estimates in the multilevel model. Due to missing observations when an individual does not appear in a bootstrap, \code{nobs < m x B}, where \code{m < M} is the number of individuals in a given cluster, \code{M} is the total number of individuals and \code{B} is the total number of bootstrap in \code{regressboot}.}
  \item{pooled.ame}{A numeric value indicating the pooled AME, which is the mean change in cluster membership probability for a change in the level of the covariate of interest over all bootstraps and all individuals belonging to the reference cluster in the original typology.}
  \item{standard.error}{Standard error of the pooled AME, which diminishes asymptotically as the number of bootstrap increases.}
  \item{bootstrap.deviation}{The estimate for the standard deviation of the bootstrap random effect. This can be used to construct a prediction interval for the association of interest (see Roth et al. 2024 for details on how to compute this).}
  \item{individual.deviation}{The estimate for the standard deviation of the bootstrap random effect.}
  \item{bootstrap.ranef}{A vector of size \code{B} containing the estimated random effects for each bootstrap.}
  \item{individual.ranef}{A vector of size \code{m} containing the estimated random effects for each individual in the reference cluster.}
}
\references{
Roth, L., Studer, M., Zuercher, E., and Peytremann-Bridevaux, I. (2024). Robustness assessment of regressions using cluster analysis typologies: a bootstrap procedure with application in state sequence analysis. BMC medical research methodology, 24(1), 303.

Studer, M. (2013). WeightedCluster library manual: A practical guide to creating typologies of trajectories in the social sciences with R. University of Geneva.

Fernandez-Castilla, B., Maes, M., Declercq, L., Jamshidi, L., Beretvas, S. N., Onghena, P., and Van den Noortgate, W. (2019). A demonstration and evaluation of the use of cross-classified random-effects models for meta-analysis. Behavior research methods, 51(3), 1286–1304.
}
\author{
Leonard Roth
}
\note{
Requires the following packages: dplyr, DescTools, lme4
}
\seealso{
\code{\link{regressboot}}, \code{\link{rarcattables}}
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
result <- rarcat(bootout,  clustqual$clustering$cluster2, 2, "funempyes")
round(result$pooled.ame, 4)
round(result$standard.error, 4)
round(result$bootstrap.deviation, 4)
}
