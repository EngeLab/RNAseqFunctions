.matrixCheckingAndCoercion <- function(mat) {
  if(!is.data.frame(mat) & !is.matrix(mat)) {
    stop("matrix arg is not a matrix or a data.frame.")
  }
  if(is.data.frame(mat) & all(sapply(mat, is.numeric))) {
    message("Coercing data.frame to matrix.")
    mat <- as.matrix(mat)
  }
  if(is.data.frame(mat) & !all(sapply(mat, is.numeric))) {
    stop("matrix arg is a data frame with non-numeric columns.")
  }
  return(mat)
}

.checkNarg <- function(n, cpm) {
  if(!is.numeric(n)) {
    stop("n nust be a length 1 numeric.")
  }
  if(n > nrow(cpm)) {
    warning(paste0("n > nrow(matrix), setting n to ", nrow(cpm), "."))
    n <- min(n, nrow(cpm))
  }
  return(n)
}

.has_zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

.check0range <- function(vec) {
  if(.has_zero_range(vec)) {
    warning(paste0("All rows have an identical measured metric ", vec[1], "."))
  }
}
