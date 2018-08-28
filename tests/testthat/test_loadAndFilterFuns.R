context("countsFuns")

library(tibble)
library(dplyr)

test_that("check moveGenesToRownames", {

  #setup expected data
  e.rownames <- c("ACTB", letters[1:11], "ERCC-1", "__alignment_not_unique")

  #run function
  output <- moveGenesToRownames(testingCounts)

  #test
  expect_identical(rownames(output), e.rownames)
  expect_identical(output$HGN, NULL)
  expect_error(moveGenesToRownames(1:10))
  expect_error(moveGenesToRownames(matrix(1:10)))
  expect_error(moveGenesToRownames(data.frame(1:10)))
})

test_that("check convertCountsToMatrix", {

  #setup input
  input <- moveGenesToRownames(testingCounts)

  #run function
  output <- convertCountsToMatrix(input)

  #test
  expect_identical(class(output), "matrix")
  expect_error(convertCountsToMatrix(1:10))
  expect_error(convertCountsToMatrix(matrix(1:10)))
  expect_error(convertCountsToMatrix(testingCounts))
})

test_that("check removeHTSEQsuffix", {

  #setup input
  input <- moveGenesToRownames(testingCounts)

  #setup expected data
  e.colnames <- LETTERS[1:11]

  #run function
  output <- removeHTSEQsuffix(input)

  #test
  expect_identical(colnames(output), e.colnames)
  expect_error(removeHTSEQsuffix(1:10))
  expect_error(removeHTSEQsuffix(matrix(1:10)))
  expect_warning(removeHTSEQsuffix(data.frame(1:10, 1:10)))
  expect_warning(removeHTSEQsuffix(data.frame(a = 1:10, b = 1:10)))

  #setup input
  input <- paste0(LETTERS, ".htseq")

  #setup expected data
  expected <- LETTERS

  #run function
  output <- removeHTSEQsuffix(input)

  #test
  expect_identical(expected, output)

  #test input checks
  d <- data.frame(A = 1:10, B = 1:10)
  expect_error(removeHTSEQsuffix(1:10))
  expect_warning(removeHTSEQsuffix(d))
  colnames(d) <- NULL
  expect_error(removeHTSEQsuffix(d))
  colnames(d) <- paste0("V", 1:ncol(d))
  expect_warning(removeHTSEQsuffix(d))
})

test_that("check detectERCCreads", {

  #setup input
  #testingCounts only has 1 ercc read and thus generates a warning
  input <- data.frame(runif(100), row.names = c(1:8, paste0("ERCC-", 1:92)))

  #setup expected data
  expected <- c(rep(FALSE, 8), rep(TRUE, 92))

  #run function
  output <- detectERCCreads(input)

  #test
  expect_identical(output, expected)
  expect_error(detectERCCreads(1:10))
  expect_error(detectERCCreads(data.frame(1:10)))
})

test_that("check detectNonGenes", {

  #setup input
  input <- moveGenesToRownames(testingCounts)

  #setup expected data
  expected <- c(rep(FALSE, 13), TRUE)

  #run function
  output <- detectNonGenes(input)

  #test
  expect_identical(output, expected)
  expect_error(detectNonGenes(1:10))
  expect_error(detectNonGenes(data.frame(1:10)))
})

test_that("check detectLowQualityGenes", {

  #setup input
  input <- moveGenesToRownames(testingCounts)
  input <- input[!detectNonGenes(input), ]
  #testingCounts only has 1 ercc read and thus generates a warning
  input <- expect_warning(input[!detectERCCreads(input), ])

  #setup expected data
  expected <- c(rep(TRUE, 11), FALSE)
  names(expected) <- c("ACTB", letters[1:11])

  #run function
  output <- detectLowQualityGenes(input, 18)

  #test
  expect_identical(output, expected)
})

test_that("check detectLowQualityCells", {

  #setup input
  input <- moveGenesToRownames(testingCounts)
  input <- removeHTSEQsuffix(input)
  input <- input[!detectNonGenes(input), ]
  #testingCounts only has 1 ercc read and thus generates a warning
  input <- expect_warning(input[!detectERCCreads(input), ])

  #setup expected data
  expected <- c(rep(TRUE, 4), FALSE, rep(TRUE, 5), FALSE)
  names(expected) <- LETTERS[1:11]

  #run function
  #coerces to matrix from data.frame and gives message
  output <- expect_message(detectLowQualityCells(input, mincount = 30))

  #test
  expect_identical(output, expected)
  expect_error(expect_message(detectLowQualityCells(input[-1, ], mincount = 30)))
  expect_error(expect_message(detectLowQualityCells(input[-1, ], mincount = 1e5)))
})

test_that("check annotatePlate", {

  #setup expected data
  names <- c(
  "s.NJB00201.B01", "s.NJB00201.A07", "s.NJB00201.F03", "s.NJB00201.F07",
  "m.NJB00204.A09", "m.NJB00204.B05"
  )

  expected <- tibble(
  sample = names,
  plate = c(rep("NJB00201", 4), rep("NJB00204", 2))
  )

  #run function
  output <- annotatePlate(testingMeta)

  #test
  expect_identical(expected, output)
})

test_that("check annotateRow", {

  #setup expected data
  names <- c(
  "s.NJB00201.B01", "s.NJB00201.A07", "s.NJB00201.F03", "s.NJB00201.F07",
  "m.NJB00204.A09", "m.NJB00204.B05"
  )

  expected <- tibble(
  sample = names,
  row = c(2L, 1L, 6L, 6L, 1L, 2L)
  )

  #run function
  output <- annotateRow(testingMeta)

  #test
  expect_identical(expected, output)
})

test_that("check annotateColumn", {

  #setup expected data
  names <- c(
  "s.NJB00201.B01", "s.NJB00201.A07", "s.NJB00201.F03", "s.NJB00201.F07",
  "m.NJB00204.A09", "m.NJB00204.B05"
  )

  expected <- tibble(
  sample = names,
  column = c(1L, 7L, 3L, 7L, 9L, 5L)
  )

  #run function
  output <- annotateColumn(testingMeta)

  #test
  expect_identical(expected, output)
})
