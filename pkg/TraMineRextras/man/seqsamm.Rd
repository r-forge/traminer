\name{seqsamm}
\alias{seqsamm}
\alias{seqsammseq}
\alias{seqsammeha}
\alias{plot.SAMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sequence Analysis Multistate Model (SAMM) procedure%%  ~~function to do ... ~~
}
\description{
Sequence Analysis Multistate Model (SAMM) procedure aims to simultaneously study the occurrence of transitions within trajectories and the (sub-)trajectory immediately following it over a pre-defined period of time. This strategy allows including time-varying covariates in the sequence analysis framework.
}
\usage{
seqsamm(seqdata, sublength, covar = NULL)
\method{plot}{SAMM}(x, type="d", ...)
seqsammseq(samm, spell)
seqsammeha(samm, spell, typology, persper = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{seqdata}{State sequence object created with the \code{\link{seqdef}} function. The whole trajectory followed by individuals.}
  \item{sublength}{Numeric. The length of the trajectory (or subsequence) following a transition to be considered.}
  \item{covar}{Optional \code{data.frame} storing covariates of interest. These covariates are added to the final data set.}
  \item{x}{A SAMM object produced by \code{seqsamm}}
  \item{samm}{A SAMM object produced by \code{seqsamm}.}
  \item{type}{the type of the plot. Available types are "d" for state distribution plots (chronograms), "f" for sequence frequency plots, "Ht" for transversal entropy plots, "i" for selected sequence index plots, "I" for whole set index plots, "ms" for plotting the sequence of modal states, "mt" for mean times plots, "pc" for parallel coordinate plots and "r" for representative sequence plots.}
  \item{spell}{Character. The (ending) spell to consider. It should be one of the state of the \code{\link{alphabet}} of the sequences. }
  \item{typology}{Factor or character. The typology of the trajectories out of the specified ending \code{spell}. It should contain one observation per observed ending spell.}
  \item{persper}{Logical. If \code{TRUE}, the data are returned in person-period format. Otherwise, only one line per observed spell is returned.}
  \item{...}{additional plot parameters passed to \code{\link{seqplot}}.}
}
\details{
The Sequence Analysis Multistate Model (SAMM) procedure works in three steps. First, the sub-sequences over a given time span \code{sublength} following any transition in the trajectories (\code{seqdata}) are extracted. This step is achieved using the \code{seqsamm} function. Second, these sub-sequences are clustered using SA to identify typical sub-sequences of medium-term changes. This is achieved separately for each ending spell (see \code{spell} arguement). The \code{seqsammseq} function can be used to retrieve the sub-trajectories following each ending spell. Third, multistates models are used to estimate the chance (or risk) to end each type of spell by distinguishing the type of trajectory that follows (and identified by sequence analysis). This allows estimating the  effect of covariates on the chances to start each kind of sub-sequence. The \code{seqsammeha} prepare the data to estimate the competing risk models for each ending spell. Then usual competing risks models can be used. 

Generally speaking, the SAMM procedure allows studying the time spent in each state as well as the patterns of medium-term changes occurring in trajectories. The example section below provides a step by step example of how to use it.
}
\value{
The \code{seqsamm} function returns a SAMM object, storing the reorganized data as person period. 
\item{id}{Numeric. The ID of the observation as the row number in the original \code{seqdata}.}
\item{time}{Numeric. The time unit of the current observation (from the beginning of the original sequence).}
\item{begin}{Numeric. The time of the beginning of the current spell (from the beginning of the original sequence).}
\item{spell.time}{Numeric. The time elapsed from the beginning of the current spell.}
\item{transition}{Logical. Whether a transition out of the current spell occurred within this time unit.}
\item{s.1 until s.\code{sublength}...}{The state sequence following the current observation starting from 1 (current state) until \code{sublength} time units after the current observation.}
\item{lastobs}{Logical. Whether this is the last observation of the current spell, censored or not. This is useful when one want only one row per individual, for instance to plot survival curves (see example).}
\item{x}{object of class \code{SAMM} as produced by \code{seqsamm}}
\item{Optional covariate list}{The covariates provided with the \code{covar} argument.}

The function \code{seqsammseq} returns a sequence object (see \code{\link{seqdef}}) of the trajectories following an ending spell.

The function \code{seqsammeha} returns a \code{data.frame} storing the person period data of a specific ending \code{spell} (see \code{spell} argument) considering the given \code{typology} as competing risks (see \code{typology} argument). Several variables are added to the SAMM objects (see above):
\item{SAMMtypology}{factor. The events ending the specified spell using \code{"None"} when not event occurs.}
\item{SAMM...}{Logical. A logical vector specifying whether the current observation end the spell with the following \code{...} type of trajectory.}
}
\references{
Studer, M., Struffolino, E., & Fasang, A. E. (2018). Estimating the Relationship between Time-varying Covariates and Trajectories: The Sequence Analysis Multistate Model Procedure. \emph{Sociological Methodology}, 48(1), 103–135. \doi{10.1177/0081175017747122}
}
\author{
Matthias Studer
}
\examples{
data(mvad)
mvad.seq <- seqdef(mvad, 17:86)

## For sake of simplicity we recode all "education" states to only one common state.
mvad.seq  <- seqrecode(mvad.seq, list("education"=c("FE", "HE", "school", "training")))
## We now have three states
seqdplot(mvad.seq)

###########################################################################
##  STEP I: Subsequence extraction
###########################################################################

## We start by extracting all subsequence of length 6
## We also add covariates from the mvad data frame
mvad.samm <- seqsamm(mvad.seq, 6, covar=mvad[, c("Grammar", "funemp", "gcse5eq")])
## Plotting the results to check what are the transitions out of each states.
plot(mvad.samm)
summary(mvad.samm)


###########################################################################
### STEP II: Typology of trajectory out of joblessness
###########################################################################
## We retrieve the subsequences following a transition out of a joblessness spell
jlseq <- seqsammseq(mvad.samm, "joblessness")


## Now we create a typology of these subsequences.

## Compute the clustering using LCS
jldist <- seqdist(jlseq, method="LCS")
## For sake of simplicity, keep only 2 groups
library(cluster)
jlclust <- pam(jldist, diss=TRUE, k=2, cluster.only=TRUE)
## Specify some names for our typology (here joblessness1 or joblessness2).
jltype <- paste0("joblessness", jlclust)


###########################################################################
### STEP III: Competing risks model of trajectories out of joblessness
###########################################################################

## Get the data to estimate competing risks models of the kind of trajectory
## out of jobjlessness
## We specify the SAMM object, the ending spell (joblessness) and our typology. 
jleha <- seqsammeha(mvad.samm, "joblessness", jltype)

\dontrun{
## Now jleha stores the data in person period format for competing risks 
## Discrete time model using multinomial regression
library(nnet)
multinom(SAMMtypology~spell.time+Grammar+funemp+gcse5eq, data=jleha)

## We can also have only one line per ending spell
## Plot the results
library(survival)
jleha <- seqsammeha(mvad.samm, "joblessness", jltype, persper=FALSE)
plot(survfit(Surv(spell.time, SAMMjoblessness1)~1, data=jleha))
## Cox model
summary(coxph(Surv(spell.time, SAMMjoblessness1)~gcse5eq+Grammar+funemp, data=jleha))
## Most of the time methods for recurrent events should be used.
## See for instance the coxme library to do so.

library(coxme)
summary(coxme(Surv(spell.time, SAMMjoblessness1)~gcse5eq+Grammar+funemp+(1|id), data=jleha))
}

###########################################################################
### Now repeat steps II and III for employment and then education
### (Not shown here)
###########################################################################

}
