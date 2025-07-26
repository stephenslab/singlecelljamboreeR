#' @title Fit flashier NMF model by Initializing with NNLM
#'
#' @description Since the \dQuote{greedy} initialization approach
#'   implemented in flashier does not always work well for
#'   non-negative matrix factorization (NMF), here we implement two
#'   other strategies that will often work better: (1) cycle between
#'   several rounds of greedy initialization followed by backfitting
#'   to achieve the desired number of factors; or (2) get a
#'   \dQuote{good} initialization using the NNLM package, then use
#'   this good initialization to fit an NMF using
#'   flashier. Sometimes both (1) and (2) will produce similar
#'   results, and in other cases the results might be quite
#'   different. The first strategy is implemented with
#'   \code{greedy_init = TRUE}; the second strategy is implemented
#'   with \code{greedy_init = FALSE}.
#' 
#' @param data The data matrix. It may be a sparse or dense matrix.
#'   However, if NNLM is used (\code{greedy_init = FALSE}), the
#'   initialization step will convert the sparse matrix to a dense
#'   matrix (and will issue a warning when doing so).
#'
#' @param k The desired number of factors in the matrix
#'   factorization.
#'
#' @param greedy_init If \code{greedy_init = TRUE}, use the
#'   \dQuote{greedy initialization} strategy, which performs several
#'   rounds of greedy initialization followed by backfitting until
#'   the desired number of factors is reached. If \code{greedy_init
#'   = FALSE}, use NNLM to search for a good initialization, then
#'   run the flashier backfitting to improve on this initialization.
#'
#' @param max_greedy_cycles The maximum number of cycles of greedy
#'   initialization to perform. This is only a safeguard to avoid
#'   the situation where it runs for forever.
#' 
#' @param maxiter The maximum number of backfitting iterations in each
#'  of the calls to \code{\link[flashier]{flash_backfit}}.
#'
#' @param n.threads The number of threads/CPUs used in
#'   \code{\link[NNLM]{nnmf}} (if it is used).
#'
#' @param verbose This is the \dQuote{verbose} setting for both
#'   \code{\link[flashier]{flash_backfit}} and, if called,
#'   \code{\link[NNLM]{nnmf}}.
#' 
#' @param \dots Additional arguments passed to
#'   \code{\link[flashier]{flash_init}}.
#' 
#' @return A \sQuote{flash} object.
#'
#' @examples
#' library(flashier)
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
#' @import Matrix
#' @importFrom stats runif
#' @importFrom ebnm ebnm_point_exponential
#' @importFrom flashier flash_init
#' @importFrom flashier flash_factors_init
#' @importFrom flashier flash_backfit
#' 
#' @export
#' 
flashier_nmf <- function (data, k, greedy_init = TRUE, 
                          max_greedy_cycles = 10, maxiter = 100, 
                          n.threads = 1, verbose = 1, ...) {

  # Get the size of the data matrix.
  n <- nrow(data)
  m <- ncol(data)

  if (greedy_init) {

    # USE GREEDY INITIALIZATION
    # -------------------------
    #
    # TO DO.
    #
  } else {

    # INITIALIZE USING NNLM
    # ---------------------
    #
    # TO DO: Check that NNLM is installed.
    # 
 
    # First, get a rough rank-1 NMF using NNLM.
    if (is.matrix(data)) {
      data_dense <- data
    } else {
      warning("Converting data to a (dense) matrix; this dense matrix ",
              "may exceed memory limits")
      data_dense <- as.matrix(data)
    }
  

    # Second, get a rough rank-k NMF using NNLM.
    init <- NNLM::nnmf(data_dense,k = 1,loss = "mse",method = "scd",
                       max.iter = 10,verbose = verbose,
                       n.threads = n.threads)
    W0 <- cbind(init$W,matrix(runif(n*(k-1)),n,k-1))
    H0 <- rbind(init$H,matrix(runif(m*(k-1)),k-1,m))
    nmf <- NNLM::nnmf(data_dense,k,init = list(W = W0,H = H0),loss = "mse",
                      method = "scd",max.iter = 10,verbose = verbose,
                      n.threads = n.threads)

    # Third, refine the rank-k NMF using flashier.
    out <- flash_init(data,...)
    out <- flash_factors_init(out,list(nmf$W,t(nmf$H)),
                              ebnm_point_exponential)
    out <- flash_backfit(out,extrapolate = FALSE,maxiter = maxiter,
                         verbose = verbose)
    out <- flash_backfit(out,extrapolate = TRUE,maxiter = maxiter,
                         verbose = verbose)
  }

  return(out)
}
