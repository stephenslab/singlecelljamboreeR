#' @title Title Goes Here
#'
#' @description Description goes here.
#'
#' @param fl Describe the fl input argument here.
#'
#' @param kset Describe the kset input argument here.
#'
#' @param dat Describe the dat input argument here.
#'
#' @param \dots Additional parameters.
#' 
#' @return Describe the return value here.
#'
#' @importFrom ebnm ebnm_point_exponential
#' @importFrom ebnm ebnm_point_laplace
#' @importFrom flashier flash_init
#' @importFrom flashier flash_factors_init
#' 
#' @export
#'
convert_factors_nn_to_pn <- function (fl, kset, dat, ...) {
  k <- fl$n_factors
  out <- flash_init(dat,...)

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
