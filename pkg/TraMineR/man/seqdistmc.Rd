\name{seqdistmc}
\alias{seqdistmc}
\title{Multichannel distances between sequences}
\description{Compute multichannel pairwise optimal matching (OM) distances between sequences by deriving the substitution costs from the costs of the single channels.
Works with OM and its following variants:
distance based on longest common subsequence (LCS),
Hamming distance (HAM), and Dynamic Hamming distance (DHD).}

\usage{
seqdistmc(channels, method, norm="none", indel=1, sm=NULL,
     with.missing=FALSE, full.matrix=TRUE, link="sum", cval=2,
     miss.cost=2, cweight=NULL, what="diss") }

\arguments{
 \item{channels}{A list of state sequence objects defined
    with the \code{\link{seqdef}} function, each state sequence object corresponding
    to a "channel".}

 \item{method}{a character string indicating the metric to be used.
    One of \code{"OM"} (Optimal Matching), \code{"LCS"} (Longest Common Subsequence),
    \code{"HAM"} (Hamming distance), \code{"DHD"} (Dynamic Hamming distance).}

  \item{norm}{
    String.
    Default: \code{"none"}.
    The normalization method to use. See \code{\link{seqdist}}.
    }

 \item{indel}{A vector with an insertion/deletion cost for each channel (OM method).}

 \item{sm}{A list with a substitution-cost matrix for each channel (OM, HAM and DHD
   method) or a list of method names for generating the substitution-costs
   (see \code{\link{seqsubm}}).}

 \item{with.missing}{Must be set to \code{TRUE} when sequences contain non deleted
   gaps (missing values) or when channels are of different length. See details.}

 \item{full.matrix}{If \code{TRUE} (default), the full distance matrix is returned.
   If \code{FALSE}, an object of class \code{\link{dist}} is returned.}

 \item{link}{One of \code{"sum"} or \code{"mean"}. Method to compute the "link"
   between channels. Default is to sum the substitution costs.}

 \item{cval}{Substitution cost for \code{"CONSTANT"} matrix, see \code{\link{seqsubm}}.}

 \item{miss.cost}{Missing values substitution cost, see \code{\link{seqsubm}}.}

 \item{cweight}{A vector of channel weights. Default is 1 (same weight for each channel).}
 
 \item{what}{Character string. What output should be returned? One of \code{"diss", "sm", "seqmc"}.}
}

\details{The \code{seqdistmc} function builds a state sequence by combining the channels, derives the multichannel indel and substitution costs from the indel and substitution costs of each channel (following the strategy proposed by \cite{Pollock, 2007}), and computes the multichannel distances using the multichannel distances. The available metrics (see 'method' option) are optimal matching (\code{"OM"}), longest common subsequence (\code{"LCS"}), Hamming distance (\code{"HAM"}), and Dynamic Hamming Distance (\code{"DHD"}). See \code{\link{seqdist}} for more information about distances between sequences.

Normalization may be useful when dealing with sequences that are not all of the same length. For details on the applied normalization, see \code{\link{seqdist}}. }

\value{When \code{what="diss"}, a matrix of pairwise distances between multichannel sequences.\cr
When \code{what="sm"}, the matrix of substitution costs with three attributes: \code{indel} the indel, \code{alphabet} the alphabet of the combined state sequences, and \code{cweight} the channel weight used.\cr
When \code{seqmc}, the combined state sequence object.
}

\seealso{ \code{\link{seqsubm}}, \code{\link{seqdef}}, \code{\link{seqdist}}. }

\references{ Pollock, Gary (2007) Holistic trajectories: a study of combined employment, housing and family careers by using multiple-sequence analysis. \emph{Journal of the Royal Statistical Society: Series A} \bold{170}, Part 1, 167--183. }

\examples{data(biofam)

## Building one channel per type of event left, children or married
bf <- as.matrix(biofam[, 10:25])
children <-  bf==4 | bf==5 | bf==6
married <- bf == 2 | bf== 3 | bf==6
left <- bf==1 | bf==3 | bf==5 | bf==6

## Building sequence objects
child.seq <- seqdef(children)
marr.seq <- seqdef(married)
left.seq <- seqdef(left)

## Using transition rates to compute substitution costs on each channel
mcdist <- seqdistmc(channels=list(child.seq, marr.seq, left.seq),
 	method="OM", sm =list("TRATE", "TRATE", "TRATE"))

## Using a weight of 2 for children channel and specifying substitution-cost
smatrix <- list()
smatrix[[1]] <- seqsubm(child.seq, method="CONSTANT")
smatrix[[2]] <- seqsubm(marr.seq, method="CONSTANT")
smatrix[[3]] <- seqsubm(left.seq, method="TRATE")
mcdist2 <- seqdistmc(channels=list(child.seq, marr.seq, left.seq),
	method="OM", sm =smatrix, cweight=c(2,1,1)) }

\author{Matthias Studer (with Gilbert Ritschard for the help page)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dissimilarity measures}
