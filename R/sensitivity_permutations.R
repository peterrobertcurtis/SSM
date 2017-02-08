# This file contains necessary functions for performing permutation tests on
# SSM.  There is a function to identify all effects and interactions.  There is
# a function to construct the relevant Q matrices for a given SSM.  There is a
# function to generate an augmented lattice design around a design.  There is a
# function to run a permutation test on a given Q matrix.

#' Construct the Q matrix for an effect or interaction in an SSM.
#'
#' This function takes an ssm as an input along with a vector of integers that
#' denotes the effect or interaction of interest.  It then computes the relevant
#' Q matrix that is used to compute the FANOVA variance of that
#' effect/interaction using the formulation \eqn{y^TQy}.
#'
#' @param ssm An SSM class object.
#' @param indices A vector of integers identifying the effect or interaction of
#'   interest.  For example, \code{indices = c(1, 2)} indicates the interaction
#'   between \eqn{x_1} and \eqn{x_2} while \code{indices = 1} indicates the main
#'   effect of \eqn{x_1}.
#' @param type (optional) This is passed to \code{\link{constructH}} and
#'   determines the method of computing the \eqn{H} matrix.  Must be 1 or 2.
#'   used to construct \eqn{Q}. See \code{\link{constructH}} for details.
#' @export
#' @return The Q matrix used in the variance formulation \eqn{y^TQy}.
constructQ <- function(ssm, indices, type = 1){
  if (length(ssm@total_variance < 1)) ssm <- update.sensitivity(ssm)
  if(length(indices) == 1){
    Delta <- diag(constructDelta(ssm) * (ssm@main_ind[, indices]))
    H <- constructH(ssm, type)
    Q <- H %*% Delta %*% t(H)
  } else {
    ind <- apply(ssm@basis[ , factors] > 0 ,1, prod)
    if(length(factors) != ssm@dimension){
      not_factors <- setdiff(1:ssm@dimension, factors)
      not_ind <- apply(as.matrix(ssm@basis[ , not_factors], ncol=length(not_factors)) == 0, 1, prod)
      ind <- ind * not_ind
    }
    Delta <- diag(constructDelta(ssm) * ind)
    H <- constructH(ssm, type)
    Q <- H %*% Delta %*% t(H)
  }
  return(Q)
}

#' Construct the H matrix defining the indicator polynomials used in the
#' construction of an SSM.
#'
#' Given an SSM this function computes the H matrix that defines a size
#' \eqn{n} basis of the model such that \eqn{Hy=\theta}.
#'
#' The method of computation when the \code{type} argument is set to 2 assumes
#' that the first \eqn{d+1} columns of the design model matrix are associated
#' with the affine polynomials in the model basis.  Therefore this method may
#' not be appropriate if a user-defined basis has been used and \code{type}
#' should be set to 1.
#'
#' @param ssm An SSM class object.
#' @param type (optional) Determines the way \eqn{H} is computed. The default
#'   \code{type = 1} will construct and invert C, taking the appropriate
#'   sub-matrix. However, \code{type = 2} will use submatrices of the design
#'   model matrix and the K matrix which may or may not be quicker (to be
#'   determined through testing).
#' @export
#' @return The H matrix.
constructH <- function(ssm, type = 1){
  if(type == 1){
    C <- rbind(cbind(diag(0, ssm@design_size), ssm@design_model_matrix),
               cbind(t(ssm@design_model_matrix), ssm@K))
    H <- solve(C)[1:ssm@design_size, -(1:ssm@design_size)]
    return(H)
  } else if (type == 2){
    X0   <- ssm@design_model_matrix[ , 1:(ssm@dimension+1)]
    X1   <- ssm@design_model_matrix[ , -(1:(ssm@dimension+1))]
    Kinv <- solve(ssm@K[-(1:(ssm@dimension+1)), -(1:(ssm@dimension+1))])
    Ainv <- solve(X1 %*% Kinv %*% t(X1))
    Binv <- solve(t(X0) %*% Ainv %*% X0)
    H1   <- Binv %*% t(X0) %*% Ainv
    H2   <- Kinv %*% t(X1) %*% (Ainv - Ainv %*% X0 %*% Binv %*% t(X0) %*% Ainv)
    H <- t(rbind(H1, H2))
    return(H)
  } else stop("The argument 'type' must be either 1 or 2.")
}

#' Construct the full \eqn{\Delta} matrix.
#'
#' This constructs the diagonal \eqn{\Delta} matrix that contains the
#' normalisation constants for the basis polynomials.
#'
#' This assumes that the basis is formed from Legendre polynomials and is not
#' appropriate if the basis is user-defined.
#'
#' @param ssm An SSM class object.
#'
#' @return A vector containing the diagonal entries of the \eqn{\Delta} matrix.
#' @export
constructDelta <- function(ssm){
  if (ssm@dimension == 1) {
    Delta <- as.numeric(2 / (2 * ssm@basis + 1))
  } else {
    Delta <- apply(2 / (2 * ssm@basis + 1), 1, prod)
  }
  return(Delta)
}
