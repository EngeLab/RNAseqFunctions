context("featureSelection")

####
#nTopVar
####

test_that("check nTopVar with expected input", {

  #setup input data
  set.seed(934)
  mat <- matrix(sample(1:100, 10), ncol = 2)

  #setup expected data
  expected <- c(5L, 2L)

  #run function
  output <- nTopVar(mat, 2)

  #test
  expect_identical(expected, output)
})

test_that("check nTopVar with bad n class", {
  #test
  expect_error(nTopVar(matrix(1:10, ncol = 2), "A"))
})

test_that("check nTopVar with n too large", {
  #test
  expect_warning(nTopVar(matrix(1:10, ncol = 2), 20))
})

test_that("check nTopVar with identical varience", {
  #test
  expect_warning(nTopVar(matrix(rep(1:10, each = 2), ncol = 2), 2))
})
