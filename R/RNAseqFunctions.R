#' A package for Enge lab RNAseq analysis.
#'
#' A package for Enge lab RNAseq analysis.
#'
#' @name RNAseqFunctions
#' @docType package
#' @author Author: Martin Enge and Jason T. Serviss
#'
#' @useDynLib RNAseqFunctions
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

