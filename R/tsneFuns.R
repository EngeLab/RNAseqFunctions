
#' pearsonsCor
#'
#' Calculates 1-Pearson's correlation between columns and returns results as a
#' dist object.
#'
#' @name pearsonsCor
#' @rdname pearsonsCor
#' @author Jason T. Serviss
#' @param counts The matrix holding expression values.
#' @param select Optional; row indexes to include.
#' @keywords pearsonsCor
#' @examples
#' pearsonsCor(testingCounts[, -1])
#'
NULL

#' @export
#' @importFrom stats as.dist cor

pearsonsCor <- function(cpm, select = NULL) {
  if(is.null(select)) select <- 1:nrow(cpm)
  cpm <- .matrixCheckingAndCoercion(cpm)
  as.dist(1 - cor(cpm[select, ], method = "p"))
}

#' runTsne
#'
#' This method runs t-SNE based on an input of 1-Pearson's correlation. Uses
#' the \code{\link[Rtsne]{Rtsne}} function from the Rtsne package.
#'
#' @name runTsne
#' @rdname runTsne
#' @aliases runTsne
#' @param my.dist dist; Typically produced with the \code{\link{pearsonsCor}}
#'  function.
#' @param dims integer; Argument to \code{\link[Rtsne]{Rtsne}}.
#'  Output dimensionality (default: 2)
#' @param theta numeric; Argument to \code{\link[Rtsne]{Rtsne}}. Speed/accuracy
#'  trade-off (increase for less accuracy), set to 0.0 for exact TSNE
#'  (default: 0)
#' @param initial_dims integer; Argument to \code{\link[Rtsne]{Rtsne}}. The
#'  number of dimensions that should be retained in the initial PCA step
#'  (default: 50)
#' @param max_iter integer; Argument to \code{\link[Rtsne]{Rtsne}}. Number of
#'  iterations (default: 2000)
#' @param perplexity numeric; Argument to \code{\link[Rtsne]{Rtsne}}. Perplexity
#'  parameter (should not be bigger than 3 * perplexity < nrow(X) - 1).
#'  (default = 10)
#' @param seed The desired seed to set before running.
#' @param is_distance logical; Argument to \code{\link[Rtsne]{Rtsne}}. Indicate
#'  whether X is a distance matrix (default = TRUE).
#' @param ... Additional arguments to pass on
#' @return Matrix containing the new representations for the objects.
#' @author Jason T. Serviss
#' @keywords runTsne
#' @examples
#' pc <- pearsonsCor(testingCounts[, -1])
#' t <- runTsne(pc, perplexity = 2)
#'
NULL

#' @rdname runTsne
#' @importFrom Rtsne Rtsne
#' @export

runTsne <- function(
  my.dist, dims = 2, theta = 0, initial_dims = 50, max_iter = 2000,
  perplexity = 10, seed = 11, is_distance = TRUE, ...
){
  if(!is.matrix(my.dist) & class(my.dist) != "dist") {
    stop("The my.dist arg is not a martix or a dist object.")
  }
  if(class(my.dist) == "dist" & !is_distance) {
    warning("my.dist is class dist and is_distance is FALSE")
  }
  if(is.data.frame(my.dist)) {
    my.dist <- .matrixCheckingAndCoercion(my.dist)
  }

  set.seed(seed)
  my.tsne <- Rtsne(
    my.dist, dims = dims, initial_dims = initial_dims, max_iter = max_iter,
    perplexity = perplexity, theta = theta, is_distance = is_distance
  )$Y
  rownames(my.tsne) <- attr(my.dist, "Labels")
  return(my.tsne)
}
