#' @title Add title here.
#'
#' @description Add description here.
#'
#' @param n Describe this argument here.
#'
#' @param verbose Describe this argument here.
#' 
#' @return Describe return value here.
#' 
#' @export
#' 
hello_world <- function (n, verbose = TRUE) {
  for (i in 1:n)
    cat("Hello world!\n")
  return(rep("Hello world!",n))
}
