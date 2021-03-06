% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GP_SSM.R
\name{concentrated.likelihood}
\alias{concentrated.likelihood}
\title{Compute the concentrated likelihood of a covariance matrix.}
\usage{
concentrated.likelihood(r, ssm)
}
\arguments{
\item{r}{A number. The length parameter of the correlation function.}

\item{ssm}{An SSM object.}
}
\value{
A number. The concentrated likelihood.
}
\description{
Computes the concentrated likelihood of the covariance matrix of an SSM
object, given a length parameter and the SSM Leave-One-Out errors.
}
\keyword{internal}

