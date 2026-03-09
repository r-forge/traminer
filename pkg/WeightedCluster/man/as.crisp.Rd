\name{as.crisp}
\alias{as.crisp}


\title{
	Conversion from Fuzzy to Crisp Clustering 
}
\description{
The function converts fuzzy clusterings into crisp ones. 
}

\usage{
as.crisp(clara)
}

\arguments{
  \item{clara}{The results of \code{seqclararange}}
}

\details{
A fuzzy clustering is converted to a crisp one by assinging each observation to the cluster showing the highest membership probability.
}
\value{
See \code{seqclararange} for details
}

\examples{

\dontrun{
	# load illustrative data
	data(mvad)

	# Creating a state sequence object
	mvad.seq <- seqdef(mvad[1:200, 17:86])

	# Fuzzy clustering

	mvadFuzzy <- seqclararange(mvad.seq, R = 3, sample.size = 20, kvals = 2:3, 
	  method="fuzzy", seqdist.args = list(method = "HAM"), 
		parallel=FALSE)
	 
	# Conversion to crisp

	mvadCrisp <- as.crisp(mvadFuzzy)

	# Plotting sequences in 3 groups

	seqdplot(mvad.seq, 
			 group = mvadCrisp$clustering$cluster3, 
			 border = NA)
}
}
