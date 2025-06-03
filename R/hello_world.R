#' @title Print Hello World
#'
#' @description Prints a \dQuote{Hello World} message to the console
#'   and provides the same message in the output.
#'
#' @param n The number of times to say \dQuote{Hello World!}
#'
#' @param verbose If \code{verbose = TRUE}, print the messages to the
#'     console.
#' 
#' @return A character vector containing the message.
#' 
#' @export
#' 
hello_world <- function (n, verbose = TRUE) {
  if (verbose) {
    for (i in 1:n)
      cat("Hello world!\n")
  }
  return(rep("Hello world!",n))
}
