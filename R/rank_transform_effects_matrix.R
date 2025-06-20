#' @title Rank-Transform Effects Matrix
#'
#' @description For n x m effects matrix X, this function transforms
#'   each column of X so that the values lie between 0 and 1, with 0
#'   being the lowest rank, and 1 being the highest. By
#'   "highest rank", we mean that the entry has the largest (most
#'   positive) value amongst all the entries in the column.
#' 
#' @param effects_matrix The n x m matrix X.
#'
#' @param compare_cols Explain here what \code{compare_cols} = TRUE does.
#'
#' @return The n x m rank-transformed effects matrix.
#' 
#' @export
#'
rank_transform_effects_matrix <- function (effects_matrix,
                                           compare_cols = FALSE) {
  effects_matrix <- apply(effects_matrix,2,rank_random_tie)
  if (compare_cols) {
    out <- effects_matrix
    m <- ncol(effects_matrix)
    for (i in 1:m) {
      x <- effects_matrix[,i]
      p <- x - apply(effects_matrix[,-i,drop = FALSE],1,max)  
      n <- x - apply(effects_matrix[,-i,drop = FALSE],1,min)
      out[,i] <- ifelse(abs(p) < abs(n),p,n)
    }
  } else {
    out <- effects_matrix
  }
  return(out)
}

# This is a helper function used by rank_transform_effects_matrix.
rank_random_tie <- function (x) {
  n <- length(x)
  return((rank(x,ties.method = "random",na.last = "keep") - 1)/(n-1))
}

