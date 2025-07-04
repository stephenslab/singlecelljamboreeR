#' @title Add Title Here
#'
#' @description Add description here.
#' 
#' @param gene_signals rows = genes, columns = signals
#'
#' @param gene_sets rows = genes, columns = gene sets
#'
#' @param gene_set_info Describe the gene_set_info input here.
#'  
#' @param min_size Describe the min_size input here.
#'
#' @param max_size Describe the max_size input here.
#'
#' @param top_genes Describe the top_genes argument here.
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
