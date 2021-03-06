% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sensitivity_SSM.R
\name{sensitivity.plot}
\alias{sensitivity.plot}
\title{Plot the sensitivity indices of a smooth supersaturated model.}
\usage{
sensitivity.plot(ssm, type = "main_total", ...)
}
\arguments{
\item{ssm}{An SSM object.  Must have relevant sensitivity indices
in the appropriate slots, either from setting \code{SA = TRUE} in
\code{\link{fit.ssm}} or by using \code{\link{update.sensitivity}}.}

\item{type}{(optional) Character. Determines the type of barplot. One of
\code{"sobol"}, \code{"main_sobol"}, \code{"main_total"}, or
\code{"total"}. Default behaviour is \code{"main_total"}.}

\item{...}{arguments passed to the \code{barplot} function call.}
}
\description{
\code{sensitivity.plot} visualises the sensitivity indices of a given smooth
supersaturated model using \code{barplot}.  If the \code{SA} flag was not set
to TRUE when \code{\link{fit.ssm}} was run to fit the model, then
\code{\link{update.sensitivity}} should be used to compute the model
variances.  If not, this function will return an error message.
}
\details{
There are four classes of plot available:
\itemize{
  \item{\code{"sobol"}}{ Produces a barplot of all Sobol indices. If there
  are more than 10 factors then Sobol indices will not have been computed for
  interactions and only the Sobol indices for main effects will be plotted.
  Main effects are in red, interactions are in grey.}
  \item{\code{"main_sobol"}}{ Produces a barplot of Sobol indices for main
  effects only.}
  \item{\code{"main_total"}}{ Produces a barplot of Total indices for main
  effects only. This is the default behaviour.}
  \item{\code{"total"}}{ Produces a barplot of Total indices for main effects
  and all second order interactions. Main effects are in red, interactions
  are in grey.}
}

Note that variables and interactions are not labelled in the plots since
there can be too many bars to label clearly.
}
\examples{
# A 20 point design in four variables
X <- matrix(runif(80, -1, 1), ncol = 4)
Y <- runif(20)
s <- fit.ssm(X, Y, SA = TRUE)
sensitivity.plot(s)

# In the next plots, the grey bars indicate interactions.
sensitivity.plot(s, "sobol")
sensitivity.plot(s, "total")
# Identifying particular indices is best done using the information held in
# the SSM object.  The following orders s@total_int_factors so the
# interaction indicated by the top row is the most important.
ind <- order(s@total_int, decreasing = TRUE)
s@total_int_factors[ind, ]
}

