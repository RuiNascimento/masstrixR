# load required libraries ======================================================
library(masstrixR)
library(MSnbase)
library(tidyverse)

# get test data ================================================================
mgf_file <- system.file("testdata",
                        "test_utilities_MS2_AutoMSn_CelegansLipids.mgf",
                        package = "masstrixR")

db_file <- system.file("testdata",
                       "test_utilities_MS2_wormLipidDb_MS2.sqlite",
                       package = "masstrixR")

# read the MS2 data ============================================================
ms2Spectra <- filterMsLevel(readAnnotatedMgfData(mgf_file), msLevel = 2)


# search for different fragmentation characteristics ===========================
# check for water loss ---------------------------------------------------------
mcols(ms2Spectra)$waterLoss <- unlist(lapply(ms2Spectra,
                                      containsNeutralLossIon,
                                      neutralLoss = 18.010565))

# check for C17 sphingoid base fragment (m/z 250.2529) -------------------------
mcols(ms2Spectra)$c17sphingoid <- unlist(lapply(ms2Spectra,
                                         containsProductIon,
                                         productIonMz = c(250.2529,
                                                          268.2635,
                                                          238.2530),
                                         multiplePi = "all"))

# perform filtering ============================================================
# filter spectra for putative sphingolipids ------------------------------------
ms2Spectra_filtered <- ms2Spectra[which(mcols(ms2Spectra)$waterLoss,
                                        mcols(ms2Spectra)$c17sphingoid)]

# filter for m/z above 200 to remove noise
ms2Spectra_filtered <- ms2Spectra_filtered[which(precursorMz(ms2Spectra_filtered) > 200)]

# perform annotation with spectra from library =================================
# adducts to search for
adducts <- c("[M+H]+", "[M+Na]+")

# empty data frame for results -------------------------------------------------
result <- tibble()

# iterate over all spectra
for(i in seq_along(ms2Spectra_filtered)) {

  # create empty Spectra object for results
  searchResult <- new("Spectra")

  # perform precursor search
  # iterate over adducts
  for(adduct in adducts) {
    searchResultClipboard <-
      searchByPrecursor(
        precursorMz(ms2Spectra_filtered[[i]]),
        db_file,
        precursorType = adduct,
        mzTol = 0.005
      )

    # if 1 or more results were found, add to searchResult
    if(length(searchResultClipboard) > 0) {
      searchResult <- append(searchResult, searchResultClipboard)
    }
  }

  result <- bind_rows(result,
                      createResultsSet(ms2Spectra_filtered[[i]],
                                       searchResult,
                                       prefix = paste0("spectrum no ", i),
                                       savePlot = FALSE))
}


# convert to tibble and print ==================================================
result <- as_tibble(result)
