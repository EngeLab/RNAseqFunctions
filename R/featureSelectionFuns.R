#' nTopVar
#'
#' Facilitates gene selection using variance.
#'
#' Returns the index for the n genes (rows) with the maximum
#' variance in the counts object. Counts per million (CPM) should be used for
#' the calculation.
#'
#' @name nTopVar
#' @rdname nTopVar
#' @aliases nTopVar
#' @param cpm matrix; Matrix containing cpm values.
#' @param n Number of genes to select.
#' @return A numeric vector containing the indices of selected genes.
#' @author Jason T. Serviss
#' @keywords nTopVar
#'
NULL

#' @rdname nTopVar
#' @importFrom stats var
#' @importFrom matrixStats rowVars
#' @export

nTopVar <- function(cpm, n) {
  n <- min(n, dim(cpm)[1])
  rv = rowVars(cpm)
  order(rv, decreasing = TRUE)[1:n]
}

#' nTopMax
#'
#' Facilitates gene selection prior to unsupervised clustering.
#'
#' Returns the index for the n genes (rows) with the maximum
#' expression in the spCounts object. The expression matrix in
#' the counts.cpm slot is used for the calculation.
#'
#' @name nTopMax
#' @rdname nTopMax
#' @aliases nTopMax
#' @param cpm matrix; Matrix containing cpm values.
#' @param n Number of genes to select.
#' @return A numeric vector containing the indices of selected genes.
#' @author Jason T. Serviss
#' @keywords nTopMax
#'
NULL

#' @rdname nTopMax
#' @export

nTopMax <- function(cpm, n) {
  n <- min(n, dim(cpm)[1])
  rv <- apply(cpm, 1, max)
  select <- order(rv, decreasing = TRUE)[1:n]
  return(select)
}

#' nTopDeltaCV
#'
#' Selected informative genes by fitting a support-vector regression to the
#' coefficient of variation (CV) as a function of the mean, and selecting genes
#' having the greatest offset from the fitted curve; this would correspond to
#' genes with higher-than-expected variance. Adpoted from Linnarsson lab python
#' implementation:
#' https://github.com/linnarsson-lab/cytograph/blob/master/cytograph/feature_selection.py
#'
#' @name nTopDeltaCV
#' @rdname nTopDeltaCV
#' @author Jason T. Serviss
#' @param cpm matrix; Matrix containing cpm values.
#' @param n Number of genes to select.
#' @keywords nTopDeltaCV
#'
#'
#' @export
NULL

#' @importFrom e1071 svm
#' @importFrom matrixStats rowMeans2 rowSds
#' @importFrom stats predict

nTopDeltaCV <- function(counts, n) {
  valid <- matrixStats::rowSums2(counts) > 0
  mu <- matrixStats::rowMeans2(counts)
  sd <- matrixStats::rowSds(counts)
  ok <- mu > 0 & sd > 0
  cv <- sd[ok] / mu[ok]

  log2_m <- log2(mu[ok])
  log2_cv <- log2(cv)

  svr_gamma <- 1000 / length(mu[ok])
  modelsvm <- svm(log2_cv ~ log2_m, gamma = svr_gamma)
  score <- log2_cv - predict(modelsvm, log2_m)
  score <- score * valid[ok]
  names(score) <- rownames(counts)[ok]
  sort(score, decreasing = TRUE)[1:n]
}
