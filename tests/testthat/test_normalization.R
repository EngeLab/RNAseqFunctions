context("normalization")

test_that("check cpm with expected input", {

  #setup normal input data
  mat <- matrix(rep(c(1e6, 1e7), each = 5), ncol = 2)

  #setup expected data
  expected <- t(t(mat) / colSums(mat) * 10^6)

  #run function
  output <- cpm(mat)

  #test
  expect_identical(expected, output)
})

test_that("check log2cpm with expected input", {

  #setup normal input data
  mat <- matrix(rep(c(1e6, 1e7), each = 5), ncol = 2)

  #setup expected data
  expected <- log2(t(t(mat) / colSums(mat) * 10^6) + 1)

  #run function
  output <- log2cpm(mat)

  #test
  expect_identical(expected, output)
})

test_that("check normalizeVec with expected input", {

  #setup normal input data
  vec <- 1:10

  #setup expected data
  expected <- c(
    0, 0.111111111111111, 0.222222222222222, 0.333333333333333,
    0.444444444444444, 0.555555555555556, 0.666666666666667, 0.777777777777778,
    0.888888888888889, 1
  )

  #run function
  output <- normalizeVec(vec)

  #test
  expect_equal(expected, output)
})

test_that("check normalizeVec with unexpected input", {
  expect_error(normalizeVec(letters))
})
