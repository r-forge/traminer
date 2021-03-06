\name{seqgen}
\alias{seqgen}
\title{Random sequences generation}
\description{
Generates random sequences.
}
\details{
Each sequence is generated by choosing a set of random numbers (with min=1 and max=length of the alphabet) using the \code{\link{runif}} function. When the probability distribution is not specified, the uniform probability distribution giving same probability to each state is used to generate the sequences.
}
\usage{
 seqgen(n, length, alphabet, p)
}
\arguments{
  \item{n}{number of sequences to generate}
  \item{length}{sequences length}
  \item{alphabet}{the alphabet from which the sequences are generated}
  \item{p}{an optional vector of probabilities for the states in the alphabet. Must be of the same length as the alphabet. If not specified, equal probabilities are used.}
}
\value{a sequence object.}
%\seealso{}
\examples{
seq <- seqgen(1000,10,1:4,c(0.2,0.1,0.3,0.4))
seqstatd(seqdef(seq))
}

\author{Alexis Gabadinho (with Gilbert Ritschard for the help page)}
\keyword{Data handling}
\keyword{State sequences} 
