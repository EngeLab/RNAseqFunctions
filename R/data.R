#' Test counts
#'
#' A small example data frame used for demonstrating and testing the counts data
#' processing functions available in the package.
#'
#' @title Simulated data used for testing and demos.
#' @docType data
#' @name testingCounts
#' @format Matrix counts, with:
#' \describe{
#' \item{rownames}{Gene names}
#' \item{colnames}{Samples/cells}
#' }
#' @keywords datasets
#' @examples
#' data(testingCounts)
"testingCounts"

#' testingMeta
#'
#' Data for testing metadata processing functions.
#'
#' @title Data for testing metadata processing functions.
#' @docType data
#' @name testingMeta
#' @format Tibble with:
#' \describe{
#' \item{sample}{Sample name}
#' \item{row}{numeric; indicates the plate row number}
#' \item{column}{numeric; indicates the plate column number}
#' }
#' @keywords datasets
#' @examples
#' data(testingMeta)
#'
"testingMeta"

#' pro.counts
#'
#' A small example dataset of procesed counts.
#'
#' @title scRNAseq data from 3 cell lines.
#' @docType data
#' @name pro.counts
#' @format Matrix counts, with:
#' \describe{
#' \item{rownames}{Gene names}
#' \item{colnames}{Samples/cells}
#' }
#' @keywords datasets
#' @examples
#' data(pro.counts)
"pro.counts"

#' pro.meta
#'
#' Metadata for the counts dataset.
#'
#' @title Metadata for the counts dataset.
#' @docType data
#' @name pro.meta
#' @format Tibble with:
#' \describe{
#' \item{sample}{Sample name}
#' \item{cellTypes}{character; indicates the cell type}
#' }
#' @keywords datasets
#' @examples
#' data(pro.meta)
#'
"pro.meta"
