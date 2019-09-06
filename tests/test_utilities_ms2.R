library(masstrixR)
library(tidyverse)
library(MSnbase)

context("MS2 utility functions")

## tests for exact mass to ion mass --------------------------------------------
test_that("Utility functions for MS2 level are correct", {

  # check that spectra are of correct length
  expect_equal(length(ms2Spectra), 10361)

  # check for correct filtering
  expect_equal(length(ms2Spectra_filtered), 166)

})
