context("metadata_general")

library(tibble)
library(dplyr)
data(testingMeta)

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
