#' @title Rank Effects Matrix
#'
#' @description Rank the effects in each column of the input matrix,
#'   from largest (rank of 1) to smallest (rank of n, where n is the
#'   number of rows).
#'
#' @param effects_matrix The input effects matrix in which rows
#'   correspond to the samples/observations and columns correspond
#'   to thedifferent signals.
#'
#' @return A matrix of the same dimension as the input matrix
#'   containing the rankings.
#'
#' @examples
#' X <- matrix(c(1:5,10:6,-(1:5)),5,3)
#' rownames(X) <- paste0("row",1:5)
#' colnames(X) <- paste0("col",1:3)
#' Y <- rank_effects(abs(X))
#' 
#' @export
#'
rank_effects <- function (effects_matrix) {
  if (!(is.matrix(effects_matrix) & is.numeric(effects_matrix)))
    stop("Input \"effects_matrix\" should be a numeric matrix")
  return(apply(-effects_matrix,2,rank))
}
