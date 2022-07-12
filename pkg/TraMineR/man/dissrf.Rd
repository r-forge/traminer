\name{dissrf}
\alias{dissrf}
\alias{seqrf}
\alias{summary.dissrf}
\alias{summary.seqrf}
\alias{print.dissrf}
\alias{print.seqrf}
\alias{plot.seqrf}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relative Frequency Groups.
}
\description{
Relative Frequency (RF) groups are equally sized groups obtained by partitioning sorted cases into \eqn{k} consecutive groups. Function \code{dissrf} returns medoid index of the RF groups and related statistics and \code{seqrf} returns in addition the RF medoid sequences.
}
\usage{
dissrf(diss,
       k=NULL,
       sortv=NULL,
       weights=NULL,
       grp.meth = "first",
       squared = FALSE)
       
seqrf(seqdata, diss,
       k=NULL,
       sortv=NULL,
       weights=NULL,
       weighted=TRUE,
       grp.meth = "first",
       squared = FALSE)
       
\method{plot}{seqrf}(x, space=0, border=NA, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{diss}{Matrix or distance object. Pairwise dissimilarities between analyzed cases.}
  \item{seqdata}{State sequence \code{stslist} object as produced by \code{\link{seqdef}}.}
  \item{k}{Integer: Number of groupings (frequency groups)}
  \item{sortv}{Real vector (of length \code{nrow(diss)}), character string, or \code{NULL}. Sorting variable or method used to compute the frequency groups. If a string, one of \code{"from.start"} and \code{"from.end"}.If \code{NULL}, the first MDS factor of \code{diss} is used. Ties are randomly ordered.}
  \item{weights}{Vector  (of length \code{nrow(diss)}) of non-negative weights. If \code{NULL} (default), equal weights are used.}
  \item{weighted}{Logical. Should weights be used when there are weights in \code{seqdata}? (default is \code{TRUE})}
  \item{grp.meth}{Character string. One of \code{"first"} (default), \code{"random"}, and \code{"prop"}. Grouping method. See details.}
  \item{squared}{Logical. Should dissimilarities be squared?}
  \item{x}{a \code{seqrf} object}
  \item{space}{space between horizontal bars representing medoid sequences (see \code{\link{barplot}})}
  \item{border}{the color to be used for the border of the bars. Use border = NA to omit borders. (see \code{\link{barplot}})}
  \item{...}{additional arguments passed to \code{\link{plot.stslist}}}
}

\details{
Function \code{dissrf} partitions the \eqn{n} cases (rows of the \code{diss} matrix) into \eqn{k} equally sized groups. First, the cases are sorted according to the \code{sortv} variable (or method). Then the groups are built by consecutively grouping the first \eqn{n/k} cases, then the next \eqn{n/k} cases, and so on.

The \code{grp.meth} argument serves when the group size (\eqn{n/k}) is not integer. With \code{grp.meth="first"}, the integer part of \eqn{n/k} is used as group size and the size of the first groups is augmented by one unit so that the sum of the group sizes equals \eqn{n}. With \code{grp.meth="random"}, randomly selected groups have their size augmented by one unit, and with \code{grp.meth="prop"}, cases at the limit between groups are assigned proportionally to each of the two groups.

For \code{seqrf}, when \code{weights=NULL} and \code{weighted=TRUE}, \code{weights} is set as the \code{weights} attribute of \code{seqdata}.

When \code{weights} is non-null, \code{grp.meth} is forced as \code{"prop"}, because it is the only method applicable to weighted data.

The function computes indicative statistics of the resulting partition, namely a pseudo R2 and a pseudo F statistics. These statistics compare the mean distance to the group medoid with the mean distance to the overall medoid. When \code{squared} is \code{TRUE}, mean squared dissimilarities are used and otherwise the R2 and F ratios are based on mean of non-squared dissimilarities. An indicative p-value of the F statistics is computed using the F distribution. This p-value should be interpreted with caution since F is not a true F value.

}

\value{\code{dissrf} returns a list of class \code{dissrfprop} when \code{grp.meth="prop"} and of class \code{dissrfcrisp} otherwise. In both cases the list also receives class \code{"dissrf"}. The elements of the list are:
  \item{medoids}{index of the group medoids}
  \item{med.names}{names (diss colnames) of the group medoids}
  \item{wg}{working matrix used by the \code{"prop"} procedure (class \code{dissrfprop} only)}
  \item{dist.list}{list with for each successive group the distances from its elements to the group medoid}
  \item{index.list}{list with for each successive group the index of its elements}
  \item{weights.list}{list with for each successive group the weights of its elements in the group}
  \item{heights}{relative group size, which may be different when \code{grp.meth} is \code{"first"} or \code{"random"}}
  \item{kmedoid.index}{vector with for each case the index of its group medoid (class \code{dissrfcrisp} only)}
  \item{kmedoid.dist}{vector with for each case the distance to its group medoid (class \code{dissrfcrisp} only)}
  \item{mdsk}{vector of group membership (class \code{dissrfcrisp} only)}
  \item{at}{positions for the boxplots of distances to group medoids}
  \item{R2}{Pseudo R2: Mean distance to the group medoids over mean distance to the overall medoid}
  \item{Fstat}{Pseudo F statistics}
  \item{pvalue}{p-value of the pseudo F (to be used with caution since F is not a true F value)}
  \item{sizes}{\code{ncase} (number of cases), \code{wsum} (sum of weights), \code{k} (number of groups), \code{gsize} (group size)}
      
\code{seqrf} returns a list of class \code{seqrfprop} when \code{grp.meth="prop"} and of class \code{seqrfcrisp} otherwise. In both cases the list also receives class \code{"seqrf"}. The elements of the list are:
    \item{seqtoplot}{RF medoid sequences as a state sequence \code{stslist} object}
    \item{rf}{the associated \code{dissrf} object} 

There are print and summary methods for objects of class \code{dissrf} and \code{seqrf}, and a plot method for objects of class \code{seqrf}
}

\seealso{\code{\link{seqrfplot}}, \code{\link{dissrep}}, and \code{\link{seqrep}}
}

\references{
Fasang, Anette Eva and Tim F. Liao. 2014. "Visualizing Sequences in the Social Sciences: Relative Frequency Sequence Plots." Sociological Methods & Research 43(4):643-676.
}
\author{
Gilbert Ritschard.
}

\examples{
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")

## Here, we use only 100 cases selected such that all elements
## of the alphabet be present.
## (More cases and a larger k would be necessary to get a meaningful example.)
biofam.seq <- seqdef(biofam[501:600, 10:25], labels=biofam.lab,
                    weights=biofam[501:600,"wp00tbgs"])
diss <- seqdist(biofam.seq, method="LCS")

## Using 12 groups and default MDS sorting
dissrf(diss=diss, k=12, grp.meth="random")

## Using 12 groups, weights, and default MDS sorting
w <- attr(biofam.seq,"weights")
dissrf(diss=diss, k=12, grp.meth="prop", weights=w)

## With a user specified sorting variable
## Here time spent in parental home: there are ties
parentTime <- seqistatd(biofam.seq)[, 1]
b.srf <- seqrf(biofam.seq, diss=diss, k=12, sortv=parentTime)
b.srf
summary(b.srf)
plot(b.srf)

}

