#' cpm
#'
#' Calculates counts per million (cpm) using a gene expression counts matrix
#' as input.
#'
#' @name cpm
#' @rdname cpm
#' @aliases cpm
#' @param counts matrix; a numeric matrix of counts.
#' @return A matrix of cpm values.
#' @author Jason T. Serviss
NULL

#' @export
#' @importFrom matrixStats colSums2

cpm <- function(counts) {
  if(is.data.frame(counts)) {
    warning("Coercing data.frame to matrix")
    counts <- as.matrix(counts)
  }
  t(t(counts) / matrixStats::colSums2(counts) * 10^6)
}

#' log.cpm
#'
#' Calculates counts per million (cpm) using a gene expression counts matrix
#' as input.
#'
#' @name cpm.log2
#' @rdname cpm.log2
#' @aliases cpm.log2
#' @param counts matrix; a numeric matrix of counts.
#' @return A matrix of log2 cpm values.
#' @author Jason T. Serviss
NULL

#' @export
#' @importFrom matrixStats colSums2

cpm.log2 <- function(counts) {
  if(is.data.frame(counts)) {
    warning("Coercing data.frame to matrix")
    counts <- as.matrix(counts)
  }
  log2(t(t(counts) / matrixStats::colSums2(counts) * 10^6) + 1)
}

#' normalizeVec
#'
#' Normalizes a vector (x) using x - min(x) / max(x) - min(x).
#'
#' @name normalizeVec
#' @rdname normalizeVec
#' @author Jason T. Serviss
#' @param vec numeric; A numeric vector.
#' @keywords normalizeVec
#' @examples
#'
#' normalizeVec(rnorm(100))
#'
#' @export

normalizeVec <- function(vec) {
  (vec - min(vec)) / (max(vec) - min(vec))
}
