# Perform a gene set enrichment analysis using susieR.
topics <- colnames(Y)
gsea <- vector("list",ncol(Y))
names(gsea) <- topics
for (i in topics) {
  cat("topic",i,"\n")
  out <- susie(X,Y[,i],L = 10,intercept = TRUE,standardize = FALSE,
               estimate_residual_variance = TRUE,refine = FALSE,
               compute_univariate_zscore = FALSE,verbose = TRUE,
               min_abs_corr = 0)
  gsea[[i]] <- out[c("KL","lbf","sigma2","V","elbo","sets","pip","alpha","mu")]
}


