library(masstrixR)
context("ion matching")

## tests for exact mass to ion mass --------------------------------------------
test_that("checked that different adduct ions are matched correctly", {

  # intra (within one ion mode) ------------------------------------------------
  adducts <- c("[M+H]+", "[M+Na]+", "[M+NH4]+")

  mz1 <- 181.070665
  mz2 <- 203.052609

  # positive mode, single charge
  expect_equal(match_mz_intra(mz1, mz2, adducts), "[M+H]+<->[M+Na]+")

  # inter (between two ion modes) ----------------------------------------------
  neg_adducts <- c("[M-H]-")
  pos_adducts <- c("[M+H]+", "[M+Na]+", "[M+NH4]+")

  mz_neg <- 179.055014
  mz_pos <- 203.052609

  # negative mode, multiple charge
  expect_equal(match_mz_inter(mz_pos, mz_neg, pos_adducts, neg_adducts), "[M+Na]+<->[M-H]-")

})
