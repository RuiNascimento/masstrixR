#' Calculate forward Dotproduct
#'
#'
#' @export
forwardDotProduct <- function(x, y, align = FALSE, mzTol = 0.005, mzTolType = "abs", binwidth = 1L, ...) {

  # no further continued
  .Deprecated(msg = "'forwardDotProduct' will be removed in the next version")

  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol, mzTolType = mzTolType)

    # calculate dot product
    dotproduct <- .dotproduct(alignedSpectra$intensity.top, alignedSpectra$intensity.bottom)

    # return result
    return(dotproduct)

  } else {

    # use binned spectra
    binnedSpectra <- bin_Spectra(x, y, ...)

    # calculate dot product
    dotproduct <- .dotproduct(binnedSpectra$intensity.top, binnedSpectra$intensity.bottom)

    # return result
    return(dotproduct)
  }
}


#' Calculate reverse Dotproduct
#'
#'
#' @export
reverseDotProduct <- function(x, y, align = FALSE, mzTol = 0.005, mzTolType = "abs", binwidth = 1L, ...) {

  # no further continued
  .Deprecated(msg = "'reverseDotProduct' will be removed in the next version")

  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol, mzTolType = mzTolType)

    # use only peals with are present in spec2
    alignedSpectra <- alignedSpectra[which(alignedSpectra$intensity.bottom >0), ]

    # calculate dot product
    dotproduct <- .dotproduct(alignedSpectra$intensity.top, alignedSpectra$intensity.bottom)

    # return result
    return(dotproduct)

  } else {

    # use binned spectra
    binnedSpectra <- bin_Spectra(x, y, ...)

    # use only peals with are present in spec2
    binnedSpectra <- binnedSpectra[which(binnedSpectra$intensity.bottom >0), ]

    # calculate dot product
    dotproduct <- .dotproduct(binnedSpectra$intensity.top, binnedSpectra$intensity.bottom)

    # return result
    return(dotproduct)
  }
}


#' Calculate common peaks
#'
#'
#' @export
commonPeaks <- function(x, y, align = FALSE, mzTol = 0.005, mzTolType = "abs", binwidth = 1L, ...) {

  # no further continued
  .Deprecated(msg = "'commonPeaks' will be removed in the next version")


  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol, mzTolType = mzTolType)

    commonPeaks <- nrow(alignedSpectra[which(alignedSpectra$intensity.top > 0 & alignedSpectra$intensity.bottom > 0),])

    # return number of common peaks
    return(commonPeaks)

  } else {

    # use binned spectra
    binnedSpectra <- bin_Spectra(x, y, ...)

    commonPeaks <- nrow(binnedSpectra[which(binnedSpectra$intensity.top > 0 & binnedSpectra$intensity.bottom > 0),])

    # return number of common peak
    return(commonPeaks)

  }
}

################################################################################
# NEW VERSIONS BELOW
################################################################################


#' Calculate forward Dotproduct
#'
#'
#' @export
forward_dotproduct <- function(x, y, align = FALSE, mzTol = 0.005, mzTolType = "abs", binwidth = 1L, ...) {

  if(align) {

    #use aligned spectra
    compared_spectra <- alignSpectra(x, y, mzTol = mzTol, mzTolType = mzTolType)

  } else {

    # use binned spectra
    compared_spectra <- bin_Spectra(x, y, ...)

  }

  dotproduct <- .dotproduct(compared_spectra$intensity.top, compared_spectra$intensity.bottom)


  return(dotproduct)

}


#' Calculate reverse Dotproduct
#'
#'
#' @export
reverse_dotproduct <- function(x, y, align = FALSE, mzTol = 0.005, mzTolType = "abs", binwidth = 1L, ...) {

  if(align) {

    #use aligned spectra
    compared_spectra <- alignSpectra(x, y, mzTol = mzTol, mzTolType = mzTolType)

  } else {

    # use binned spectra
    compared_spectra <- bin_Spectra(x, y, ...)

  }

  # use only peals with are present in spec2
  compared_spectra <- compared_spectra[which(compared_spectra$intensity.bottom >0), ]

  # calculate dot product
  dotproduct <- .dotproduct(compared_spectra$intensity.top, compared_spectra$intensity.bottom)

  return(dotproduct)
}


#' Calculate common peaks
#'
#'
#' @export
common_peaks <- function(x, y, align = FALSE, mzTol = 0.005, mzTolType = "abs", binwidth = 1L, ...) {

  if(align) {

    #use aligned spectra
    compared_spectra <- alignSpectra(x, y, mzTol = mzTol, mzTolType = mzTolType)

  } else {

    # use binned spectra
    compared_spectra <- bin_Spectra(x, y, ...)

  }

  commonPeaks <- nrow(compared_spectra[which(compared_spectra$intensity.top > 0 & compared_spectra$intensity.bottom > 0),])

  # return number of common peak
  return(commonPeaks)
}




#' Helper function to calculate dot product between two intensity vectors on same m/z scale (binned or aligned)
.dotproduct <- function(x, y) {
  as.vector(x %*% y) / (sqrt(sum(x*x)) * sqrt(sum(y*y)))
}
