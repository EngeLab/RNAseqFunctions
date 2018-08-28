context("inputChecks")

####
#.matrixCheckingAndCoercion
####

test_that("check .matrixCheckingAndCoercion with expected input", {

  #setup normal input data
  mat <- matrix(rep(1:10, each = 2), ncol = 2, byrow = TRUE)

  #setup expected data
  expected <- mat

  #run function
  output <- .matrixCheckingAndCoercion(mat)

  #test
  expect_identical(expected, output)
})

test_that("check .matrixCheckingAndCoercion with correct data.frame input", {

  #setup input data
  dat <- data.frame(A = 1:10, B = 1:10)

  #setup expected data
  expected <- matrix(rep(1:10, 2), ncol = 2, dimnames = list(NULL, c("A", "B")))

  #run function
  output <- expect_message(.matrixCheckingAndCoercion(dat))

  #test
  expect_identical(expected, output)
})

test_that("check .matrixCheckingAndCoercion with incorrect data.frame input", {

  #setup input data
  dat <- data.frame(A = 1:10, B = 1:10, C = LETTERS[1:10])

  #test
  expect_error(.matrixCheckingAndCoercion(dat))
})

####
#.checkNarg
####

test_that("check .checkNarg with expected input", {

  #setup normal input data
  mat <- matrix(rep(1:10, each = 2), ncol = 2, byrow = TRUE)
  n <- 2

  #setup expected data
  expected <- n

  #run function
  output <- .checkNarg(n, mat)

  #test
  expect_identical(expected, output)
})

test_that("check .checkNarg with unexpected class", {
  #test
  expect_error(.checkNarg(matrix(1:10), "A"))
})

test_that("check .checkNarg with n > nrow(cpm)", {

  #setup normal input data
  mat <- matrix(rep(1:10, each = 2), ncol = 2, byrow = TRUE)
  n <- 11

  #setup expected data
  expected <- 10

  #run function
  output <- expect_warning(.checkNarg(n, mat))

  #test
  expect_identical(expected, output)
})

####
#.has_zero_range
####

test_that("check .has_zero_range with expected input", {
  #test
  expect_identical(.has_zero_range(rep(0, 10)))
})
