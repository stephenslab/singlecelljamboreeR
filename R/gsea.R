#' @title Perform Gene-set Enrichment Analyses Using SuSiE
#'
#' @description Add description here.
#' 
#' @param gene_signals A matrix of signals in which rows correspond to
#'   genes and columns correspond to the different signals. If there
#'   is only one signal, this can be a vector.
#'
#' @param gene_sets A binary matrix (sparse or dense) of gene sets in
#'   which rows correspond to genes and columns correspond to gene
#'   sets. An entry of 1 means the gene is a member of the gene set.
#'
#' @param gene_set_info Optional data frame containing additional
#'   information about the gene sets.
#'  
#' @param min_size Gene sets with fewer genes than this will not be
#'     considered.
#'
#' @param max_size Gene sets with more genes than this will not be
#'   considered.
#'
#' @param top_genes The number of \dQuote{top genes} to include in the 
#' 
#' @param L The maximum number of selected gene sets; passed as the
#'   \dQuote{L} argument to \code{\link[susieR]{susie}}.
#'     
#' @param coverage The \dQuote{coverage} parameter in
#'   \code{\link[susieR]{susie}}.
#'
#' @param max_iter The \dQuote{max_iter} parameter in
#'   \code{\link[susieR]{susie}}.
#'
#' @param tol The \dQuote{tol} parameter in
#'   \code{\link[susieR]{susie}}.
#' 
#' @param verbose If \code{verbose = TRUE}, print updates about 
#'   progress of the analysis.
#'
#' @return A list containing (a) the SuSiE model fits (one for each
#'   signal); (b) a tibble containing the gene set enrichment results.
#'
#' @examples
#' 
#' # Add an example here illustrating the use of perform_gsea().
#' library(pathways)
#' set.seed(1)
#' data(gene_sets_human)
#' i <- which(!is.na(gene_sets_human$gene_info$Ensembl))
#' y <- gene_sets_human$gene_sets[i,"M973"] + 0.1 * rnorm(length(i))
#' out <- perform_gsea(y,gene_sets_human$gene_sets[i,])
#' out$selected_gene_sets
#'
#' @importFrom Matrix colSums
#' @importFrom susieR susie
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' 
#' @export
#' 
perform_gsea <- function (gene_signals, gene_sets, gene_set_info = NULL,
                          min_size = 10, max_size = 400, top_genes = 10, 
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

  # Verify and process the "gene_set_info" input.
  if (!is.null(gene_set_info)) {
    if (!is.data.frame(gene_set_info))
      stop("Input argument \"gene_set_info\" should be a data frame")
    if (ncol(gene_sets) != nrow(gene_set_info))
      stop("Input argument \"gene_set_info\" should have one column for each row of \"gene_sets\"")
    if (is.null(rownames(gene_set_info)))
      stop("Input argument \"gene_set_info\" should have named rows")
    if (any(colnames(gene_sets) != rownames(gene_set_info)))
      stop("The rows of \"gene_set_info\" should have the same names as the columns of \"gene_sets\"")     
  }

  # Align the gene sets with the gene statistics.
  genes        <- intersect(rownames(gene_signals),rownames(gene_sets))
  gene_signals <- gene_signals[genes,,drop = FALSE]
  gene_sets    <- gene_sets[genes,,drop = FALSE]
  i            <- match(genes,rownames(gene_signals))
  j            <- match(genes,rownames(gene_sets))
  gene_signals <- gene_signals[i,,drop = FALSE]
  gene_sets    <- gene_sets[j,,drop = FALSE]

  # Next, remove gene sets with fewer than min_size genes and with
  # more than max_size genes.
  x <- colSums(gene_sets)
  i <- which(x >= min_size & x <= max_size)
  gene_sets     <- gene_sets[,i,drop = FALSE]
  gene_set_info <- gene_set_info[i,]

  # Perform a gene set enrichment analysis using susieR.
  if (verbose) {
    cat("Number of gene signals:",ncol(gene_signals),"\n")
    cat("Number of gene sets:",ncol(gene_sets),"\n")
    cat("Number of genes:",nrow(gene_sets),"\n")
  }
  signals <- colnames(gene_signals)
  susie_fits <- vector("list",ncol(gene_signals))
  names(susie_fits) <- signals
  for (i in signals) {
    if (verbose)
      cat("signal ",i,":\n",sep = "")
    susie_fits[[i]] <- 
      susie(gene_sets,gene_signals[,i],L = L,intercept = TRUE,
            standardize = FALSE,estimate_residual_variance = TRUE,
            estimate_prior_variance = TRUE,refine = FALSE,
            estimate_prior_method = "EM",coverage = coverage,
            max_iter = max_iter, tol = tol,
            compute_univariate_zscore = FALSE,
            verbose = verbose,min_abs_corr = 0)
  }

  # Compile the enriched gene sets into a single table.
  if (verbose)
    cat("Compiling the results.\n")
  selected_gene_sets <- NULL
  for (i in signals) {
    out <- compile_gsea_table(susie_fits[[i]],gene_signals[,i],gene_sets,
                              gene_set_info,top_genes)
    out <- bind_cols(tibble(signal = rep(i,nrow(out))),out)
    selected_gene_sets <- bind_rows(selected_gene_sets,out)
  }

  # Make a couple adjustments to the final table.
  selected_gene_sets$signal <- factor(selected_gene_sets$signal)
  selected_gene_sets$CS     <- factor(selected_gene_sets$CS)
  return(list(selected_gene_sets = selected_gene_sets,
              susie_fits = susie_fits))
}

# Generate a data frame summarizing the results of the SuSiE-based
# gene set enrichment analysis.
#
#' @importFrom Matrix colSums
#' @importFrom tibble tibble
#' @importFrom tibble as_tibble
#' @importFrom dplyr bind_cols
#' 
compile_gsea_table <- function (s, gene_signal, gene_sets, gene_set_info, 
                                top_genes = 10) {

  # Initialize the output.
  out <- NULL

  # Get the credible sets.
  cs <- s$sets$cs

  # Get the gene set identifiers and the gene set sizes.
  gene_set_ids   <- colnames(gene_sets)
  gene_set_sizes <- colSums(gene_sets)

  # Add labels to some of the susie outputs.
  n                 <- length(s$lbf)
  names(s$lbf)      <- paste0("L",1:n)
  rownames(s$alpha) <- paste0("L",1:n)
  rownames(s$mu)    <- paste0("L",1:n)

  # Reorder the CSs by Bayes factor.
  css <- names(cs)
  i   <- order(s$lbf[css],decreasing = TRUE)
  css <- css[i]

  # Repeat for each CS.
  for (i in css) {

    # Get information about the variables (i.e., the gene sets)
    # included in the CS.
    j <- cs[[i]]
    n <- length(j)
    x <- tibble(CS        = rep(i,n),
                gene_set  = gene_set_ids[j],
                lbf       = rep(s$lbf[i],n),
                pip       = s$alpha[i,j],
                coef      = s$mu[i,j],
                genes     = gene_set_sizes[j],
                top_genes = vector("list",n))
    if (!is.null(gene_set_info))
      x <- bind_cols(x,as_tibble(gene_set_info[j,]))

    # For each gene set in the CS, extract the "top genes", which are
    # the genes in the gene set with the largest signals (in
    # magnitude).
    t <- 1
    for (k in j) {
      genes <- which(gene_sets[,k] > 0)
      y     <- abs(gene_signal[genes])
      genes <- genes[order(y,decreasing = TRUE)]
      if (length(genes) > top_genes)
        genes <- genes[seq(1,top_genes)]
      x$top_genes[[t]] <- names(genes)
      t <- t + 1
    }

    # Reorder the gene sets by PIP.
    rows <- order(x$pip,decreasing = TRUE)
    x    <- x[rows,]
    out  <- rbind(out,x)
  }

  rownames(out) <- NULL
  return(out)
}
