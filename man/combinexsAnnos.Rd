\name{combinexsAnnos}
\alias{combinexsAnnos}
\title{Check CAMERA ion species annotation due to matching with opposite ion mode}
\description{
 This function check annoations of ion species with the help of a sample from opposite ion mode.
 As first step it searches for pseudospectra  from the positive and the negative sample
 within a retention time window. For every result the m/z differences between both samples 
 are matched against specific rules, which are combinations from pos. and neg. ion species.
 As example M+H and M-H with a m/z difference of 2.014552. If two ions matches such a difference,
 the ion annotations are changed (previous annotation is wrong), confirmed or added.
 Returns the peaklist from one ion mode with recalculated annotations.
 }
\usage{
  combinexsAnnos(xsa.pos, xsa.neg, pos=TRUE, tol=2, ruleset=NULL)
}
\arguments{
  \item{xsa.pos}{xsAnnotate object with positive ion mode}
  \item{xsa.neg}{xsAnnotate object with neagtive ion mode}
  \item{pos}{If TRUE the peaklist from the positive mode is returned, 
  if FALSE the negative}
  \item{tol}{Retention time window in seconds}
  \item{ruleset}{Matrix of matching rules, see example}
}

\details{
  Both xsAnnotate object should be full processed (grouping and annotation). 
  Without previous annotation the resulting peaklist only includes annotation 
  with matches peaks from both mode according to the rule(s). With ruleset=NULL 
  the function only looks for M+H/M-H pairs.
  The ruleset is a two column matrix with includes rule indices from the rule table
  of both xsAnnotate objects. 
  ruleset <- cbind(1,1) would create the M+H/M-H rule, since the first rule of xsa.pos@ruleset
  and xsa.neg@ruleset is M+H respectively M-H.
  Only rules with identical charge can be combined!
}

\value{
  Returns a (normal) CAMERA peaklist with a additional column neg. Mode or pos. Mode,
  where matching peaks from the opposite mode are noted.
}

\examples{
 \dontrun{
 #Searches for M+H/M-H combinations within a retention time window of 2 seconds
 peaklist.pos <- combinexsAnnos(xsa.pos, xsa.neg, tol=2)
 }
}

\author{Carsten Kuhl <ckuhl@ipb-halle.de>}
\keyword{methods}
