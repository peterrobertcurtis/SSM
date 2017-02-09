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
#'   effect of \eqn{x_1}. If set to zero, the output matrix will define the
#'   total variance.
#' @param type (optional) This is passed to \code{\link{constructH}} and
#'   determines the method of computing the \eqn{H} matrix.  Must be 1 or 2.
#'   used to construct \eqn{Q}. See \code{\link{constructH}} for details.
#' @export
#' @return The Q matrix used in the variance formulation \eqn{y^TQy}.
constructQ <- function(ssm, indices, type = 1){
  if (length(ssm@total_variance) < 1) ssm <- update.sensitivity(ssm)
  if(indices == 0){
    Delta <- diag(constructDelta(ssm))
    Delta[1,1] <- 0
  } else if(length(indices) == 1){
    Delta <- diag(constructDelta(ssm) * (ssm@main_ind[, indices]))
  } else {
    ind <- apply(ssm@basis[ , indices] > 0 ,1, prod)
    if(length(indices) != ssm@dimension){
      not_factors <- setdiff(1:ssm@dimension, indices)
      not_ind <- apply(as.matrix(ssm@basis[ , not_factors], ncol=length(not_factors)) == 0, 1, prod)
      ind <- ind * not_ind
    }
    Delta <- diag(constructDelta(ssm) * ind)
  }
  H <- constructH(ssm, type)
  Q <- H %*% Delta %*% t(H)
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

#' Identify which main effect a permutation is associated with
#'
#' This function looks at a matrix and determines which columns are not changed
#' by a permutation. TRUE indicates that the column values are unchanged.
#'
#' @param X A matrix.
#' @param perm A vector of integers denoting a permutation of design points.
#'
#' @return A vector of integers
#' @export
checkPermutedDesign <- function(X, perm){
  as.logical(apply(X[perm, ] == X, 2, prod))
}

#' Return the levels in a column of a matrix
#'
#' This returns the unique values in a specified column of a matrix.
#'
#' This is used by \code{\link{permutableIndicesByCol}}.
#'
#' @param X A matrix.
#' @param col A number indicating the column.
#' @export
columnLevels <- function(X, col){
  unique(X[, col])
}

#' Identify the positions of a given value in a vector
#'
#' This is used by \code{\link{permutableIndicesByCol}}.
#' @param i A value to search for.
#' @param vec The vector to search.
#' @export
levelIndices <- function(i, vec){
  which(vec == i)
}

#' Identify permutable matrix rows by column.
#'
#' This returns a list where each list entry is a set of matrix rows that can
#' be permuted without changing the column denoted by \code{col}.
#' @param X A matrix.
#' @param col A number indicating the column that should remain unchanged.
#' @export
permutableIndicesByCol <- function(X, col){
  L   <- columnLevels(X, col)
  res <- lapply(L, levelIndices, vec = X[ , col])
  return(res)
}


#' Generate possible permutation blocks for a matrix column
#'
#' This generates a nested list containing all possible permutation blocks that
#' leave a specified column unchanged.
#'
#' @param X A matrix.
#' @param col A number indicating the column that is to be unchanged by
#'   permutation.
#' @export
permutationBlocksByCol <- function(X, col){
  P <- permutableIndicesByCol(X, col)
  L <- lapply(P, combinat::permn)
  return(L)
}

#' Find the number of possible permutations for each block.
#'
#' This uses the output of \code{\link{permutationBlocksByCol}} and finds the
#' number of permuatations in each permutation block.
#'
#' @param permBlocks A nested list containing the possible permutation
#'   blocks generated by \code{\link{permutationBlocksByCol}}.
#' @export
blockLengths <- function(permBlocks){
   sapply(permBlocks, length)
}

#' Generate single random integer between 1 and N
#'
#' @param N An integer.
#' @export
randomOneToN <- function(N){
  sample(1:N, 1)
}

#' Get nth permutation from a permutation block.
#'
#' @param permBlock A single permutation block.
#' @param n An integer.
#' @export
getPermutationFromNumber <- function(permBlock, n){
  permBlock[[n]]
}

#' Get the identity permutation from the i-th block
#'
#' This is used by \code{\link{createPermutationByCol}} to put the generated
#' permutation in the correct order.
#'
#' @param i A number.
#' @param permBlocks A nested list containing the possible permutation
#'   blocks generated by \code{\link{permutationBlocksByCol}}.
#' @export
origPerm <- function(i, permBlocks){
  permBlocks[[i]][[1]]
}

#' Create a permutation that only permutes inside the allowed blocks
#'
#' This function operates on the permutation blocks list created by
#' \code{\link{permutationBlocksByCol}} and the vector of block lengths created
#' by \code{\link{blockLengths}}.  It is used by
#' \code{\link{createPermutationByCol}}.
#'
#' @param permBlocks A nested list containing the possible permutation
#'   blocks generated by \code{\link{permutationBlocksByCol}}.
#' @param bl A vector usually generated by \code{\link{blockLengths}} that
#'   gives the number of permutations in each block.
#' @return A vector indicating a row permutation.
#' @export
createPermutation <- function(permblocks, bl){
  p  <- sapply(bl, randomOneToN)
  perm <- unlist(mapply(getPermutationFromNumber, permblocks, p, SIMPLIFY = FALSE))
  orig <- unlist(sapply(1:length(permblocks), origPerm, permBlocks = permblocks,
                        simplify = FALSE))
  return(perm[orig])
}

#' Create n permutations that do not change the i-th column of a matrix
#'
#' WARNING - this uses a random sampling to generate permutations so it can be
#' inefficient. Once many permutations have been generated, the chance of
#' a sample being a new permutation decreases.  As a check, the generation will
#' end if \code{exit} permutations are generated in a row without finding a new
#' one.
#'
#' @param X A matrix
#' @param i The column to remain unchanged
#' @param n The number of permutations desired.
#' @param exit (optional) If this many permutations are generated in a row
#'  without finding a new one, the function will exit and return the
#'  permutations already generated.
#'
#' @return A list of permutations
#' @export
createPermutationByCol <- function(X, i, n, exit = 1000){
  P  <- permutationBlocksByCol(X, i)
  bl <- blockLengths(P)
  res <- vector("list", n)
  p <- NULL
  for (j in 1:n){
    count <- 1
    while (any(sapply(res, identical, p)) && count < exit){
      p <- createPermutation(P, bl)
      count <- count + 1
      if (count == exit){
        print(paste("No new permutation found after", count,
                    "tries. Ending generation."))
        return(res[1:(j-1)])
      }
    }
    res[[j]] <- p
  }
  return(res)
}

