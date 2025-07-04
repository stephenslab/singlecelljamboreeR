#' @title Add Title Here
#'
#' @description Add description here.
#' 
#' @param gene_signals rows = genes, columns = signals
#'
#' @param gene_sets rows = genes, columns = gene sets
#'
#' @param min_size Describe the min_size input here.
#'
#' @param max_size Describe the max_size input here.
#'
#' @param L Describe the L argument here.
#'
#' @param coverage Describe the coverage argument here.
#'
#' @param max_iter Describe the max_iter argument here.
#'
#' @param tol Describe tol argument here.
#' 
#' @param verbose Describe the verbose argument here.
#'
#' @return Describe the return value here.
#'
#' @examples
#' 
#' # Add an example here illustrating the use of perform_gsea().
#'
#' @importFrom Matrix colSums
#' @importFrom susieR susie
#' 
#' @export
#' 
perform_gsea <- function (gene_signals, gene_sets, 
                          min_size = 10, max_size = 400,
                          L = 10, coverage = 0.95, max_iter = 100,
                          tol = 0.001, verbose = TRUE) {

  # Verify and process the "gene_signals" input.
  gene_signals <- as.matrix(gene_signals)
  k <- ncol(gene_signals)
  if (is.null(rownames(gene_signals)))
    stop("Input argument \"gene_signals\" should be a matrix with ",
         "named rows")
  if (is.null(colnames(gene_signals)))
    colnames(gene_signals) <- paste0("k",1:k)

  # Verify and process the "gene_sets" input.
  if (!((is.matrix(gene_sets) & is.numeric(gene_sets)) | 
        inherits(gene_sets,"dgCMatrix")))
    stop("Input argument \"gene_sets\" should be a matrix ",
         "(or a sparse matrix)")
  if (is.null(rownames(gene_sets)) |
      is.null(colnames(gene_sets)))
    stop("Input argument \"gene_sets\" should have named rows and ",
         "named columns")

  # Align the gene sets with the gene statistics.
  genes        <- intersect(rownames(gene_signals),rownames(gene_sets))
  gene_signals <- gene_signals[genes,,drop = FALSE]
  gene_sets    <- gene_sets[genes,,drop = FALSE]
  i <- match(genes,rownames(gene_signals))
  j <- match(genes,rownames(gene_sets))
  gene_signals <- gene_signals[i,,drop = FALSE]
  gene_sets    <- gene_sets[j,,drop = FALSE]

  # Next, remove gene sets with fewer than min_size genes and with
  # more than max_size genes.
  x <- colSums(gene_sets)
  i <- which(x >= min_size & x <= max_size)
  gene_sets <- gene_sets[,i,drop = FALSE]

  # Perform a gene set enrichment analysis using susieR.
  signals <- colnames(gene_signals)
  gsea <- vector("list",ncol(gene_signals))
  names(gsea) <- signals
  for (i in signals) {
    if (verbose)
      cat("signal ",i,":\n",sep = "")
    gsea[[i]] <- susie(gene_sets,gene_signals[,i],L = L,intercept = TRUE,
                       standardize = FALSE,estimate_residual_variance = TRUE,
                       estimate_prior_variance = TRUE,refine = FALSE,
                       estimate_prior_method = "EM",coverage = coverage,
                       max_iter = max_iter, tol = tol,
                       compute_univariate_zscore = FALSE,
                       verbose = verbose,min_abs_corr = 0)
  }

  return(list(gsea = gsea,gene_sets = gene_sets))
}
