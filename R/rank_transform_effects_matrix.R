#' @title Rank-Transform Effects Matrix
#'
#' @description With \code{compare_cols = FALSE}, transforms each
#'   column of X so that the values lie between 0 and 1, with 0
#'   being the lowest rank, and 1 being the highest. By
#'   \dQuote{highest rank}, we mean that the entry has the largest
#'   (most positive) value amongst all the entries in the
#'   column. With \code{compare_cols = TRUE}, the \dQuote{least
#'   extreme} relative ranks are obtained by comparing to the ranks
#'   of the other columns.  This results in values between -1 (when
#'   the rank is 0 in the current column and 1 in all other columns)
#'   and +1 (when the rank is 1 in the current column and 0 in all
#'   other columns). When the current column has identical rank to
#'   one or more columns, the \dQuote{least extreme} relative rank
#'   is zero.
#' 
#' @param effects_matrix The n x m matrix X.
#'
#' @param compare_cols When \code{compare_cols = TRUE}, instead of
#'   returning the ranks, the \dQuote{least extreme} relative ranks
#'   are returned.
#'
#' @param compare_dims Describe compare_dims argument here.
#'
#' @return The n x m rank-transformed effects matrix.
#' 
#' @export
#'
rank_transform_effects_matrix <-
  function (effects_matrix,
            compare_cols = FALSE,
            compare_dims = seq(1,ncol(effects_matrix))) {
  # effects_matrix <- apply(effects_matrix,2,rank_random_tie)
  if (compare_cols) {
    out <- effects_matrix
    for (i in colnames(effects_matrix)) {
      x <- effects_matrix[,i]
      j <- setdiff(compare_dims,i)
      out[,i] <- apply(x - effects_matrix[,j],1,
                       function (x) x[which.min(abs(x))])
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

