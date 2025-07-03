#' @title Add Title Here
#'
#' @description Add description here.
#' 
#' @param gene_signals Describe the gene_signals input here.
#'
#' @param gene_sets
#'
#' @return Describe the return value here.
#'
#' @examples
#' # Add an example here illustrating the use of
#' # perform_gsea().
#' 
#' @export
#' 
perform_gsea <- function (gene_signals, gene_sets, 
                          min_size = 10, max_size = 400) {

  # Verify and process the "gene_signals" input.
  gene_signals <- as.matrix(gene_signals)
  k <- ncol(gene_signals)
  if (is.null(rownames(gene_signals)))
    stop("Input argument \"gene_signals\" should be a matrix with ",
         "named rows")
  if (is.null(colnames(gene_signals)))
    colnames(gene_signals) <- paste0("k",1:k)

  # Verify and process the "gene_sets" nput.
  if (!((is.matrix(x) & is.numeric(x)) | inherits(x,"dgCMatrix")))
    stop("Input argument \"gene_sets\" should be a matrix ",
         "(or a sparse matrix)")
  if (is.null(rownames(gene_sets)) |
      is.null(colnames(gene_sets)))
    stop("Input argument \"gene_sets\" should have named rows and ",
         "named columns")

  # Align the gene-set data with the gene-wise statistics.
genes <- intersect(rownames(Y),gene_sets_mouse$gene_info$Symbol)
Y <- Y[genes,]
i <- match(genes,gene_sets_mouse$gene_info$Symbol)
X <- X[i,]
rownames(X) <- rownames(Y)
gene_sets_mouse$gene_info <- gene_sets_mouse$gene_info[i,]

# Next, remove gene sets with fewer than 10 genes and with more than
# 400 genes. Gene sets with a large number of genes are less likely to
# be interesting, and slow down the enrichment analysis, so they are
# removed.
i <- which(colSums(X) >= 10 & colSums(X) <= 400)
X <- X[,i]


  return(NULL)
}
