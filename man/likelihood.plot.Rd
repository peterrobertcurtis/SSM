% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GP_SSM.R
\name{likelihood.plot}
\alias{likelihood.plot}
\title{Plot the concentrated likelihood of an SSM.}
\usage{
likelihood.plot(ssm, xrange = c(0, 1000), grid = 200)
}
\arguments{
\item{ssm}{An SSM object.}

\item{xrange}{(optional) The range of the x axis. Set to \code{c(0, 1000)}
by default.}

\item{grid}{(optional) A number. The number of points used to plot the curve.}
}
\description{
Plot the concentrated likelihood used to estimate the parameters of the
metamodel error estimating Gaussian process.
}
\details{
As a diagnostic it can be helpful to look at the concentrated likelihood
function of the correlation function used to estimate the metamodel error.
Flat likelihood functions make it difficult to pick a suitable \code{r}
length parameter. Note that \code{r} and \code{sigma} can be set manually.
}
\examples{
data(attitude)
X <- transform11(attitude[ 2:7])
Y <- attitude[ , 1]
s <- fit.ssm(X, Y, GP = TRUE)
likelihood.plot(s)
likelihood.plot(s, xrange = c(0, 20))
}

