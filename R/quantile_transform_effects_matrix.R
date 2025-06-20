#' @title Rank-Transform Effects Matrix
#'
#' @description For n x m effects matrix X, this function transforms
#'   each column of X so that the values lie between 0 and 1, with 0
#'   being the highest rank, and 1 being the lowest. By
#'   "highest rank", we mean that the entry has the largest (most
#'   positive) value amongst all the entries in the column (which is
#'   the opposite of the base rank function).
#' 
#' @param effects_matrix The n x m matrix X.
#'
#' @return The n x m rank-transformed effects matrix 
#' 
#' @export
#'
rank_transform_effects_matrix <- function (effects_matrix)
  apply(effects_matrix,2,rank_random_tie)

# This is a helper function used by quantile_transform_effects_matrix.
rank_random_tie <- function (x) {
  n <- length(x)
  return((rank(-x,ties.method = "random",na.last = "keep") - 1)/(n-1))
}
