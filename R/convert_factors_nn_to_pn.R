#' @title Convert Some Non-Negative Factors to Positive/Negative Factors
#'
#' @description Assuming the input flash model has non-negative
#'   factors with point-exponential priors, this function creates a
#'   new flash model in which selected factors have
#'   point-exponential priors for the columns of L and point-Laplace
#'   priors for the corresponding columns of F, which allows for the
#'   parameter estimates to be either positive or negative. Note
#'   that L and F in the new flash model are initialized to the
#'   estimates from the input flash model, so initially they are all
#'   non-negative; a subsequent call to
#'   \code{\link[flashier]{flash_backfit}} is needed to actually
#'   produce negative values.
#'
#' @param fl A flashier fit with non-negative factors, such as an
#'   output from a call to \link{flashier_nmf}.
#'
#' @param kset The indices of the non-negative factors to convert to
#'   positive/negative factors. 
#'
#' @param data The data matrix passed as an input to
#'   \code{\link[flashier]{flash_init}}.
#'
#' @param \dots Additional arguments passed to
#'   \code{\link[flashier]{flash_init}}.
#' 
#' @return A \sQuote{flash} object.
#'
#' @importFrom ebnm ebnm_point_exponential
#' @importFrom ebnm ebnm_point_laplace
#' @importFrom flashier flash_init
#' @importFrom flashier flash_factors_init
#' 
#' @export
#'
convert_factors_nn_to_pn <- function (fl, kset, data, ...) {
  out <- flash_init(data,...)
  k <- fl$n_factors
  for (i in 1:k) {
    l <- fl$L_pm[,i,drop = FALSE]
    f <- fl$F_pm[,i,drop = FALSE]
    if (is.element(i,kset)) {

      # Add a positive/negative factor. 
      out <- flash_factors_init(out,list(l,f),
                                ebnm_fn = c(ebnm_point_exponential,
                                            ebnm_point_laplace))
    } else {
        
      # Add back the non-negative factor.
      out <- flash_factors_init(out,list(l,f),
                                ebnm_fn = ebnm_point_exponential)
    }
  }

  return(out)
}
