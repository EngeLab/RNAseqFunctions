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

cpm <- function(counts) {
  counts <- .matrixCheckingAndCoercion(counts)
  cpm <- eigen_cpm(counts)
  rownames(cpm) <- rownames(counts)
  colnames(cpm) <- colnames(counts)
  cpm
}

#' log2cpm
#'
#' Calculates counts per million (cpm) using a gene expression counts matrix
#' as input.
#'
#' @name log2cpm
#' @rdname log2cpm
#' @aliases log2cpm
#' @param counts matrix; a numeric matrix of counts.
#' @return A matrix of log2 cpm values.
#' @author Jason T. Serviss
NULL

#' @export

log2cpm <- function(counts) {
  counts <- .matrixCheckingAndCoercion(counts)
  logged <- eigen_log2cpm(counts)
  rownames(logged) <- rownames(counts)
  colnames(logged) <- colnames(counts)
  logged
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
  if(!is.numeric(vec)) stop("The input vector is not numeric.")
  (vec - min(vec)) / (max(vec) - min(vec))
}
