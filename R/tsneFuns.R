
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
#'
#'
#' @export
NULL

pearsonsCor <- function(cpm, select = NULL) {
  if(is.null(select)) select <- 1:nrow(cpm)
  as.dist(1 - cor(cpm[select, ], method = "p"))
}

#' runTsne
#'
#'
#' Calculates the x and y coordinates of the mean of each classified group.
#'
#'
#' This method runs t-SNE based on an input of 1-Pearson's correlation.
#'
#' @name runTsne
#' @rdname runTsne
#' @aliases runTsne
#' @param my.dist A distance object typically produced with
#'  the \code{\link{pearsonsCor}} function.
#' @param dims Argument to Rtsne. Numeric indicating the output dimensions.
#' @param theta Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param initial_dims Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param max_iter Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param perplexity Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param seed The desired seed to set before running.
#' @param is_distance Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param ... Additional arguments to pass on
#' @return A matrix containing the mean value for each gene for each
#'    classification group.
#' @author Jason T. Serviss
#' @keywords runTsne
#' @examples
#' #
#'
#'
NULL

#' @rdname runTsne
#' @importFrom Rtsne Rtsne
#' @export

runTsne <- function(
  my.dist,
  dims = 2,
  theta = 0,
  initial_dims = 50,
  max_iter = 2000,
  perplexity = 10,
  seed = 11,
  is_distance = TRUE,
  ...
){
  set.seed(seed)
  
  my.tsne <- Rtsne(
    my.dist, dims = dims, initial_dims = initial_dims, max_iter = max_iter,
    perplexity = perplexity, theta = theta, is_distance = is_distance
  )$Y
  
  rownames(my.tsne) <- attr(my.dist, "Labels")
  return(my.tsne)
}
