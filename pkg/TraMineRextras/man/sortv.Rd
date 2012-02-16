\encoding{latin1}
\name{sortv}
\alias{sortv}
\alias{sorti}

\author{Gilbert Ritschard}

\title{Sort sequences by states at the successive positions}


\description{
Returns a sorting vector which to sort state sequences in a TraMineR sequence object (\code{\link[TraMineR]{seqdef}}) by the states at the successive positions.}

\usage{
sorti(seqdata, start = "end", sort.index=TRUE)

sortv(seqdata, start = "end")

}

\arguments{
  \item{seqdata}{A state sequence object as returned by \code{\link[TraMineR]{seqdef}}.}

  \item{start}{Where to start the sort. One of \code{"beg"} (beginning) or \code{"end"}.}
  
  \item{sort.index}{Should the function return sort indexes? If \code{FALSE}, sort values are returned.}
}

\details{
    With \code{start = "end"} (default), the primary sort key is the final state, then the previous one and so on. With \code{start = "beg"}, the primary sort key is the state at the first position, then at the next one and so on.
    
    With \code{sort.index = FALSE}, the function returns a vector of values whose order will determine the wanted order. This should be used as \code{sortv} argument of the \code{\link{seqiplot}} function. With \code{sort.index = TRUE}, the function returns a vector of indexes to be used for indexing.
    
    The \code{sortv} form is an alias for \code{sorti(..., sort.index = FALSE)}.  
}

\value{
  If \code{sort.index = FALSE}, the vector of sorting values.\cr 
  Otherwise the vector of sorting indexes.
}

\seealso{
Details about \code{type = "i"} or \code{type = "I"} in
\code{\link[TraMineR]{seqplot}} .
}

\examples{
data(actcal)
actcal.seq <- seqdef(actcal[,13:24])
par(mfrow=c(1,2))
seqIplot(actcal.seq, sortv=sortv(actcal.seq), withlegend = FALSE)
seqIplot(actcal.seq, sortv=sortv(actcal.seq, start="beg"), withlegend = FALSE)
actcal.seq[sorti(actcal.seq)[1:10],]

data(mvad)
mvad.seq <- seqdef(mvad[,17:86])
par(mfrow=c(1,2))
seqIplot(mvad.seq, sortv=sortv(mvad.seq, start="end"), withlegend = FALSE)
seqIplot(mvad.seq, sortv=sortv(mvad.seq, start="beg"), withlegend = FALSE)
print(mvad.seq[sorti(mvad.seq,start="beg")[700:712],], format="SPS")

}
\keyword{util}