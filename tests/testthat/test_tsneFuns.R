context("tsneFuns")

test_that("check pearsonsCor with expected input", {

  #setup normal input data
  mat <- matrix(rep(100:110, 10), ncol = 10)

  #setup expected data
  expected <- as.dist(matrix(rep(0, 100), ncol = 10))
  attr(expected, "call") <- NULL

  #run function
  output <- pearsonsCor(mat)
  attr(output, "call") <- NULL

  #test
  expect_equal(expected, output)

  #setup normal input data
  mat <- matrix(rep(100:110, 10), ncol = 10)

  #setup expected data
  expected <- as.dist(matrix(rep(2.220446e-16, 100), ncol = 10))
  attr(expected, "call") <- NULL

  #run function
  output <- pearsonsCor(mat, 1:10)
  attr(output, "call") <- NULL

  #test
  expect_equal(expected, output)
})

test_that("check runTsne with expected input", {
  pc <- pearsonsCor(testingCounts[, -1])
  expect_silent(runTsne(pc, perplexity = 2))
  expect_error(runTsne(testingCounts[, -1], is_distance = TRUE, perplexity = 2))
  expect_warning(runTsne(pc, is_distance = FALSE, perplexity = 2))
  expect_error(runTsne(1:10))
  expect_message(runTsne(testingCounts[, -1], is_distance = FALSE, perplexity = 2))
})
