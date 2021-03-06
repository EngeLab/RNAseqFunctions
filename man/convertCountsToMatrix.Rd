% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadAndFilterFuns.R
\name{convertCountsToMatrix}
\alias{convertCountsToMatrix}
\title{convertCountsToMatrix}
\arguments{
\item{counts}{data.frame; A data frame with counts data.}
}
\value{
The counts data.frame coerced into a matrix.
}
\description{
Coerces the counts data.frame into a matrix.
}
\examples{

counts <- data.frame(a = runif(26, 1, 100), b = runif(26, 1, 100))
convertCountsToMatrix(counts)

}
\author{
Jason Serviss
}
