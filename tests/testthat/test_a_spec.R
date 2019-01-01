library(testthat)
library(civecm)

my_spec <- ci_spec()

test_that('Default specification', {
  
  expect_s3_class(my_spec, "ci_spec")
  
})