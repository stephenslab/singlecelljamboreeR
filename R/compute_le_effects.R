#' @title Compute \dQuote{Least Extreme} Effects
#'
#' @description This function is useful for identifying
#'   \dQuote{distinctive} effects; that is, effects that are very
#'   different (much larger or much smaller) than the effects in the
#'   other signals.
#'
#' @param effects_matrix The input effects matrix in which rows
#'   correspond to the samples/observations and columns correspond
#'   to thedifferent signals.
#'
#' @return A matrix of the same dimension as the input matrix
#'   containing the least extreme (l.e.) effects.
#'
#' @examples
#' X <- rbind(c(4,0,0),
#'            c(1,-1,0),
#'            c(4,0,3),
#'            c(2,0,-4),
#'            c(0,-4,0),
#'            c(0,-4,-3),
#'            c(0,-2,-4),
#'            c(1,-2,4),
#'            c(0,0,0))
#' rownames(X) <- paste0("row",1:9)
#' colnames(X) <- paste0("col",1:3)
#' Y <- compute_le_effects(X)
#' 
#' @export
#' 
compute_le_effects <- function (effects_matrix) {
  if (!(is.matrix(effects_matrix) & is.numeric(effects_matrix)))
    stop("Input \"effects_matrix\" should be a numeric matrix")
  k <- ncol(effects_matrix)
  if (k <= 1)
    return(effects_matrix)
  out <- effects_matrix

  # Repeat for each column of the effects matrix.
  for (j in 1:k) {
    x <- effects_matrix[,j]
    i <- which(x > 0)
    if (length(i) > 0) {
      y <- cbind(0,effects_matrix[i,-j])
      out[i,j] <- pmax(0,x[i] - apply(y,1,max))
    }
    i <- which(x < 0)
    if (length(i) > 0) {
      y <- cbind(0,effects_matrix[i,-j])
      out[i,j] <- pmin(0,x[i] - apply(y,1,min))
    }
  }

  return(out)
}
