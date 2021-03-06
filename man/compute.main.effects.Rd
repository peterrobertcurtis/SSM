% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sensitivity_SSM.R
\name{compute.main.effects}
\alias{compute.main.effects}
\title{Compute main effects}
\usage{
compute.main.effects(ssm)
}
\arguments{
\item{ssm}{An SSM object.}
}
\value{
An SSM object.  This is the \code{ssm} input with the
relevant slots updated.
}
\description{
This function computes the Sobol indices of the model variables, also known
as the main effects.  It identifies the relevant rows of the SSM \code{basis}
slot and sums the term variances. Called by \code{\link{update.sensitivity}}.
}
\keyword{internal}

