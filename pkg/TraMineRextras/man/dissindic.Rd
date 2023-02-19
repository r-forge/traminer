\name{dissindic}
\alias{dissindic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the marginality and gain indicators of a relationship between a covariate and a distance matrix.
}
\description{
The gain and marginality indicators allows situating cases according to a quantitative relationship between a covariate and a distance matrix. The marginality indicator quantifies the typicality of cases within each group of the explanatory covariate using a measure of distance between cases and gravity centers. The gain indicator aims to identify cases that are either illustrative of, or discordant with, the quantitative association. It is computed as the gain of information when taking the covariate into account.
}
\usage{
dissindic(diss, group, gower = FALSE, squared = FALSE, weights = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{diss}{A dissimilarity matrix or a \code{dist} object.}
  \item{group}{A categorical variable.}
  \item{gower}{Logical: Is the dissimilarity matrix already a Gower matrix?}
  \item{squared}{Logical: Should we square the provided dissimilarities?}
  \item{weights}{Optional numerical vector of case weights.}
}
\details{

The two indicators are computed within the discrepancy analysis framework (see \code{\link{dissmfacw}}). The marginality is computed as a "residual" of the discrepancy analysis. More precisely, a high value means that a sequence (or another object) is far from the center of gravity of its group, i.e. the most typical situation. A low value indicates a sequence (or another object) is close to this gravity center. 

By combining the "residuals" of the null model (without covariate) and the marginality, we can identify sequences that are better represented when using the covariate than without it. These values can also be linked to the contributions of a sequence to the between sums of squares, a concept directly linked to the explained discrepancy. The gain therefore measures the statistical gain of information for each case when taking the covariate into account.

}
\value{A \code{data.frame} containing three columns: 

	\item{group}{The categorical variable used.}
	\item{marginality}{The value of the marginality indicator.}
	\item{gain}{The value of the gain indicator .}

}
\references{
Le Roux, G., M. Studer, A. Bringé, C. Bonvalet (2023). Selecting Qualitative Cases Using Sequence Analysis: A Mixed-Method for In-Depth Understanding of Life Course Trajectories, \emph{Advances in Life Course Research}, \doi{10.1016/j.alcr.2023.100530}.

Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2011). Discrepancy analysis of state sequences, \emph{Sociological Methods and Research}, Vol. 40(3), 471-510, \doi{10.1177/0049124111415372}.

}
\author{Matthias Studer}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dissvar}} to compute a pseudo variance from dissimilarities and for a basic introduction to concepts of discrepancy analysis.
   \cr
   \code{\link{dissassoc}} to test association between objects represented by their dissimilarities and a covariate.
\cr
\code{\link{dissmfacw}} to perform multi-factor analysis of variance from pairwise dissimilarities.

}
\examples{
## Defining a state sequence object
data(mvad)
mvad <- mvad[1:100, ] ## Use only a subsample to avoid long computation time
mvad.seq <- seqdef(mvad[, 17:86])

## Building dissimilarities (any dissimilarity measure can be used)
mvad.ham <- seqdist(mvad.seq, method="HAM")

## Study association with 
di <- dissindic(mvad.ham, group=mvad$gcse5eq)

## Plot sequences sorted by gain, illustrative trajectories at the top (and counterexample at the bottom)
seqIplot(mvad.seq, group=mvad$gcse5eq, sortv=di$gain)

## Plot sequences sorted by marginality, central trajectories at the bottom.
seqIplot(mvad.seq, group=mvad$gcse5eq, sortv=di$marginality)

##Scatterplot of the indicators separated by group value (as in Le Roux, et al. 2023)
par(mfrow=c(1, 2))
## Plot for the "no" category
plot(di$gain[mvad$gcse5eq=="no"], di$marginality[mvad$gcse5eq=="no"], main="No gcseq5q", 
	xlim=range(di$gain), ylim=range(di$marginality))
abline(h=mean(di$marginality), v=0) ## Draw reference lines
plot(di$gain[mvad$gcse5eq=="yes"], di$marginality[mvad$gcse5eq=="yes"], main="Yes gcseq5q", 
	xlim=range(di$gain), ylim=range(di$marginality))
abline(h=mean(di$marginality), v=0) ## Draw reference lines

}

