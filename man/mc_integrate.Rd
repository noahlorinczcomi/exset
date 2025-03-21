% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{mc_integrate}
\alias{mc_integrate}
\title{Monte Carlo Integration}
\usage{
mc_integrate(lam, niter = 10000, alpha = 0.05)
}
\arguments{
\item{lam}{eigenvalues of positive-definite LD matrix}

\item{niter}{number of Monte Carlo replicates used to approximate integral}

\item{alpha}{type I error rate}
}
\description{
This internal function allows you to integrate over a half plane in the joint space of independent weighted chi-squares
}
\examples{
mc_integrate()
}
\keyword{integration}
