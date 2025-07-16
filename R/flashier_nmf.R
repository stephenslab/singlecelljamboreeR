#' @title Title Goes Here.
#'
#' @description Add description of the function here.
#'
#' @param data Describe the data input here.
#'
#' @param k Describe the k input here.
#'
#' @param maxiter Describe the maxiter input here.
#'
#' @param n.threads Describe the n.threads input here.
#'
#' @param verbose Describe the verbose input here.
#' 
#' @param \dots Additional arguments passed to ...
#' 
#' @return Describe the return value here.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' m <- 200
#' k <- 3
#' L <- matrix(runif(n*k),n,k)
#' F <- matrix(runif(m*k),m,k)
#' E <- matrix(runif(n*m,0,0.01),n,m)
#' X <- tcrossprod(L,F) + E
#' fl <- flashier_nmf(X,k = 3)
#' L_est <- ldf(fl,type = "i")$L
#' F_est <- ldf(fl,type = "i")$F
#' ks <- c(1,3,2)
#' plot(L,L_est[,ks],pch = 20,xlab = "true",ylab = "estimated")
#' abline(a = 0,b = 1,lty = "dotted",col = "magenta")
#' plot(F,F_est[,ks],pch = 20,xlab = "true",ylab = "estimated")
#' abline(a = 0,b = 1,lty = "dotted",col = "magenta")
#'
#' @importFrom stats runif
#' @importFrom NNLM nnmf
#' @importFrom ebnm ebnm_point_exponential
#' @importFrom flashier flash_init
#' @importFrom flashier flash_factors_init
#' @importFrom flashier flash_backfit
#' 
#' @export
#' 
flashier_nmf <- function (data, k, maxiter = 100, n.threads = 1, 
                          verbose = 1, ...) {

  # Get the size of the data matrix.
  n <- nrow(data)
  m <- ncol(data)

  # First, get a rough rank-1 NMF using NNLM.
  init <- nnmf(data,k = 1,loss = "mse",method = "scd",
               max.iter = 10,verbose = verbose,
               n.threads = n.threads)

  # Second, get a rough rank-k NMF using NNLM.
  W0 <- cbind(init$W,matrix(runif(n*(k-1)),n,k-1))
  H0 <- rbind(init$H,matrix(runif(m*(k-1)),k-1,m))
  nmf <- nnmf(data,k,init = list(W = W0,H = H0),loss = "mse",
              method = "scd",max.iter = 10,verbose = verbose,
              n.threads = n.threads)

  # Third, refine the rank-k NMF using flashier.
  out <- flash_init(data,...)
  out <- flash_factors_init(out,list(nmf$W,t(nmf$H)),ebnm_point_exponential)
  out <- flash_backfit(out,extrapolate = FALSE,maxiter = maxiter,
                       verbose = verbose)
  out <- flash_backfit(out,extrapolate = TRUE,maxiter = maxiter,
                       verbose = verbose)
  return(out)
}
