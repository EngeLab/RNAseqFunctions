#' moveGenesToRownames
#'
#' Moves the first column of the counts data.frame to rownames and removes
#' the old column.
#'
#' @name moveGenesToRownames
#' @rdname moveGenesToRownames
#' @param counts data.frame; A data frame with counts data.
#' @return The counts data.frame is returned with the first column moved to
#' the rownames of the data.frame.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(LETTERS, a = runif(26, 1, 100))
#' moveGenesToRownames(counts)
#'
NULL
#' @export

moveGenesToRownames <- function(counts) {
  if(!length(dim(counts))) {
    stop("dim(counts) must have a positive length.")
  }
  if(class(counts) != "data.frame") {
    stop("Counts is not a data.frame")
  }
  if((class(counts[[1]]) != "character") & (class(counts[[1]]) != "factor")) {
    stop("The first column of counts is not class character or class factor.")
  }
  rownames(counts) <- counts[[1]]
  counts[[1]] <- NULL
  return(counts)
}

#' convertCountsToMatrix
#'
#' Coerces the counts data.frame into a matrix.
#'
#' @name convertCountsToMatrix
#' @rdname convertCountsToMatrix
#' @param counts data.frame; A data frame with counts data.
#' @return The counts data.frame coerced into a matrix.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(a = runif(26, 1, 100), b = runif(26, 1, 100))
#' convertCountsToMatrix(counts)
#'
NULL
#' @export

convertCountsToMatrix <- function(counts) {
  if(class(counts) != "data.frame") {
    stop("Counts is not a data.frame")
  }
  if(!all(sapply(counts, class) %in% c("numeric", "integer"))) {
    stop("Non-numeric columns detected. Should gene names be moved to rownames?")
  }
  as.matrix(counts)
}

#' removeHTSEQsuffix
#'
#' HTSeq adds the suffix ".htseq" to column names when it reports counts. This
#' function removes that suffix from the column names of the supplied counts
#' data.frame.
#'
#' @name removeHTSEQsuffix
#' @rdname removeHTSEQsuffix
#' @param data data.frame or character; A data frame with sample IDs as colnames
#'  or a vector of sample IDs.
#' @return The counts data.frame with the ".htseq" suffix removed from the
#'  column names.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(a.htseq = runif(26), b.htseq = runif(26))
#' removeHTSEQsuffix(counts)
#'
NULL
#' @export

removeHTSEQsuffix <- function(data) {
  .inputChecks_removeHTSEQsuffix(data)
  if(class(data) == "data.frame") {
    .dataframeChecks_removeHTSEQsuffix(data)
    colnames(data) <- gsub("(.*)\\.htseq$", "\\1", colnames(data))
    return(data)
  }
  if(class(data) == "character") {
    data <- gsub("(.*)\\.htseq$", "\\1", data)
    return(data)
  }
}

.inputChecks_removeHTSEQsuffix <- function(data) {
  if(!class(data) %in% c("data.frame", "character")) {
    stop("The data argument must be a character vector or data.frame.")
  }
}

.dataframeChecks_removeHTSEQsuffix <- function(data) {
  if(is.null(colnames(data))) {
    stop("is.null(colnames(data)) returned TRUE.")
  }
  if(all(colnames(data) == paste0("V", 1:ncol(data)))) {
    warning("Your colnames are V1..Vi. Are these sample names?" )
  }
  if(!any(grepl("\\.htseq", colnames(data)))) {
    warning("Could not find the .htseq suffix in colnames(data)")
  }
}

#' detectERCCreads
#'
#' Detects the 92 ERCC reads in the rownames of the counts data.frame. ERCC
#' reads must be named with the standard naming convention that matches the
#' regex "^ERCC\\-[0-9]*$".
#'
#' @name detectERCCreads
#' @rdname detectERCCreads
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @param regex Regular expression to extract ERCC reads.
#' @param warn Logical indicating if a warning should be issued if < 92 ERCC
#'  reads are detected.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row contains an ERCC read. A warning is issued if all 92
#' ERCC reads are not detected.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(runif(100), row.names = c(1:8, paste0("ERCC-", 1:92)))
#' detectERCCreads(counts)
#'
NULL
#' @export

detectERCCreads <- function(counts, regex = "^ERCC\\-[0-9]*$", warn = TRUE) {
  if(class(rownames(counts)) != "character") {
    stop("rownames(counts) is not of class character.")
  }
  if(all(rownames(counts) == as.character(1:nrow(counts)))) {
    m <- "rownames(counts) = 1:nrow(counts). Are gene names in rownames counts?"
    stop(m)
  }
  ercc <- grepl(regex, rownames(counts))

  if(sum(ercc) != 92 & warn) {
    warning("Couldn't detect all ERCC reads.")
  }

  return(ercc)
}

#' detectNonGenes
#'
#' HTSeq outputs several "non-genes" in the counts output. These include:
#' "__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned",
#' "__alignment_not_unique". Since these are measurements of the quantification
#' quality and are often undesirable in downstream analysis, this function
#' detects them in the submitted counts data.frame to allow easy removal.
#'
#' @name detectNonGenes
#' @rdname detectNonGenes
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row contains a non-gene.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(runif(10), row.names = c(1:9, "__no_feature"))
#' detectNonGenes(counts)
#'
NULL
#' @export

detectNonGenes <- function(counts) {
  if(class(rownames(counts)) != "character") {
    stop("rownames(counts) is not of class character.")
  }
  if(all(rownames(counts) == as.character(1:nrow(counts)))) {
    m <- "rownames(counts) = 1:nrow(counts). Are gene names in rownames counts?"
    stop(m)
  }

  nonGenes <- c(
    "__no_feature", "__ambiguous", "__too_low_aQual",
    "__not_aligned", "__alignment_not_unique"
  )
  rownames(counts) %in% nonGenes
}

#' detectLowQualityGenes
#'
#' In gene expression counts data if can often be the case that some genes are
#' not detected. This can simply be due to the fact that the gene is not
#' expressed in the tissue or, in addition, that the sequencing depth was not
#' sufficient to detect the gene. In addition, some genes may be detected but in
#' so few samples or at such a low level that it makes the quantified value
#' highly unreliable. In these cases, it is desireable to remove the gene before
#' downstream analysis which is facilitated by this function.
#'
#' @name detectLowQualityGenes
#' @rdname detectLowQualityGenes
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @param mincount numeric; A minimum rowSum for which rows with a higher rowSum
#' will be detected. Default = 0.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row meets both tested conditions.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(c(0, runif(10)), c(0, runif(10)), c(0, runif(10)))
#' detectLowQualityGenes(counts)
#'
NULL
#' @export

detectLowQualityGenes <- function(
  counts,
  mincount = 0
){
  #input checks

  bool <- rowSums(counts) > mincount
  message <- paste0(
    "Detected ", sum(!bool), " low quality genes out of ", nrow(counts),
    " genes input (", round(100 * (sum(!bool) / nrow(counts)), digits = 2),
    "%)."
  )
  print(message)
  return(bool)
}

#' detectLowQualityCells
#'
#' It is often the case that some samples from sequencing experiments are of
#' low quality, in many cases due to issues during the sample preperation stage.
#' Due to the fact that these samples represent a high level of technical noise,
#' it is often desirable to remove these before downstream analysis which is
#' facilitated by this function. The function achieves this using two methods.
#' First, the mincount argument detects samples whose sum across all genes is >
#' mincount. Second, we utilize a house keeping gene and assume its expression
#' to be normally distributed. We then detect samples where the probability of
#' the expression for the house keeping gene in that sample is greater than the
#' quantile.cut argument.
#'
#' @name detectLowQualityCells
#' @rdname detectLowQualityCells
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames and sample names as colnames.
#' @param mincount numeric; A minimum colSum for which columns with a higher
#' colSum will be detected. Default = 4e5.
#' @param geneName character; The gene name to use for the quantile cutoff. This
#' must be present in the rownames of the counts argument. Default is ACTB.
#' @param quantileCut numeric; This indicates probability at which the quantile
#' cutoff will be calculated using the normal distribution. Default = 0.01.
#' @return A logical vector with length = ncol(counts) that is TRUE when the
#' counts data.frame column contains a sample with colSums > mincount.
#' @author Jason Serviss
#' @examples
#' c <- moveGenesToRownames(testingCounts)[1:12, ]
#' detectLowQualityCells(c, geneName = "ACTB", mincount = 30)
#'
NULL
#' @export
#' @importFrom stats median qnorm

detectLowQualityCells <- function(
  counts,
  mincount = 4e5,
  geneName = 'ACTB',
  quantileCut = 0.01
){
  #input checks
  ##check counts matrix
  counts <- .matrixCheckingAndCoercion(counts)

  ##check that geneName is in rownames counts
  if(!geneName %in% rownames(counts)) {
    stop("geneName is not found in rownames(counts)")
  }

  #setup output vector
  output <- vector(mode = "logical", length = ncol(counts))
  names(output) <- colnames(counts)

  #colsums check
  colsums <- colSums(counts)
  cs <- colsums > mincount
  output[cs] <- TRUE

  if(sum(cs) < 2) {
    stop("One or less samples passed the colSums check.")
  }

  #house keeping check
  counts.log <- cpm.log2(counts[, cs])
  cl.act <- counts.log[geneName, ]
  cl.act.m <- median(cl.act)
  cl.act.sd <- sqrt(
    sum((cl.act[cl.act > cl.act.m] - cl.act.m) ^ 2) /
    (sum(cl.act > cl.act.m) - 1)
  )
  my.cut <- qnorm(p = quantileCut, mean = cl.act.m, sd = cl.act.sd)
  bool <- counts.log[geneName, ] > my.cut
  output[cs] <- cs[cs] & bool

  message <- paste0(
  "Detected ", sum(!output), " low quality cells out of ", ncol(counts),
  " cells input (", round(100 * (sum(!output) / ncol(counts)), digits = 2),
  "%)."
  )
  print(message)
  return(output)
}

#' Annotate plate.
#'
#' Uses the standard sample naming nomenclature to add plate row to metadata.
#' Sample names should follow: (s|m)\\.platename\\.platePosition.
#' platePosition is in the form row and column without a space where row is a
#' LETTER (A-H) and column is a number (1-12).
#'
#' @name annotatePlate
#' @rdname annotatePlate
#' @param data tibble; A tibble containing the sample names using standard
#'  nomenclature in a column named "sample".
#' @return The metadata tibble with an additional column indicating plate name.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr mutate
#' @importFrom stringr str_replace

annotatePlate <- function(data) {
  mutate(data, plate = str_replace(sample, "^^.\\.([A-Z0-9]*)\\....", "\\1"))
}

#' Annotate row.
#'
#' Uses the standard sample naming nomenclature to add plate row to metadata.
#' Sample names should follow: (s|m)\\.platename\\.platePosition.
#' platePosition is in the form row and column without a space where row is a
#' LETTER (A-H) and column is a number (1-12).
#'
#' @name annotateRow
#' @rdname annotateRow
#' @param data tibble; A tibble containing the sample names using standard
#'  nomenclature in a column named "sample".
#' @return The metadata tibble with an additional column indicating row as a
#' numeric value.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr "%>%" mutate select
#' @importFrom stringr str_replace

annotateRow <- function(data) {
  rowPos <- NULL
  data %>%
  mutate(rowPos = str_replace(sample, "^.\\.[A-Z0-9]*\\.(.)..", "\\1")) %>%
  mutate(row = match(rowPos, LETTERS[1:8])) %>%
  select(-rowPos)
}

#' Annotate column.
#'
#' Uses the standard sample naming nomenclature to add plate column to metadata.
#' Sample names should follow: (s|m)\\.platename\\.platePosition.
#' platePosition is in the form row and column without a space where row is a
#' LETTER (A-H) and column is a number (1-12).
#'
#' @name annotateColumn
#' @rdname annotateColumn
#' @param data tibble; A tibble containing the sample names using standard
#'  nomenclature in a column named "sample".
#' @return The metadata tibble with an additional column indicating column as a
#' numeric value.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr "%>%" mutate case_when select
#' @importFrom stringr str_replace

annotateColumn <- function(data) {
  colPos <- NULL
  data %>%
  mutate(colPos = str_replace(sample, "^.\\.[A-Z0-9]*\\..(..)", "\\1")) %>%
  mutate(column = case_when(
  colPos == "01" ~ 1L, colPos == "02" ~ 2L, colPos == "03" ~ 3L,
  colPos == "04" ~ 4L, colPos == "05" ~ 5L, colPos == "06" ~ 6L,
  colPos == "07" ~ 7L, colPos == "08" ~ 8L, colPos == "09" ~ 9L,
  colPos == "10" ~ 10L, colPos == "11" ~ 11L, colPos == "12" ~ 12L,
  colPos == "13" ~ 13L, colPos == "14" ~ 14L, colPos == "15" ~ 15L,
  colPos == "16" ~ 16L, colPos == "17" ~ 17L, colPos == "18" ~ 18L
  )) %>%
  select(-colPos)
}


