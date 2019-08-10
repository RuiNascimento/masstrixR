#' Checking masses from different ion modes
#'
#' This function is used to match two masses from the positive and negative ionization mode to check it they are potentially derived from the same metabolite. Based on the given adducts, all combinations are tested from positive to negative ion mode and the other way. If for both cases the mass error is below the given mass error it is assumed that the masses are derived from the same metabolite.
#'
#' @param posMz m/z value of positive ionization mode
#' @param negMz m/z value of negative ionization mode
#' @param posAdducts vector with adducts that are used for the positive ionization mode
#' @param negAdducts vector with adducts that are used for the negative ionization mode
#' @param mzTol m/z error, numeric value
#' @param mzTolType type of error used, absolute (abs) or relative (ppm)
#'
#' @examples
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{matchMassDiff}}
#' @export
match_mz_intra <- function(mz1, mz2, adducts, mzTol = 0.005, mzTolType = "abs") {

  # sanity checks
  # check if adduct definitions are good
  adduct_names <- metabolomicsUtils::get_adduct_names(mode = "all")

  if(!all(adducts %in% adduct_names) | !all(adducts %in% adduct_names)) {
    stop("one or more adducts is not fitting")
  }

  # make combinations
  adductCombinations <- as.data.frame(t(combn(adducts, 2)))

  # make df
  df <- data.frame(adduct1 = as.character(adductCombinations$V1),
                   adduct2 = as.character(adductCombinations$V2),
                   mz1 = mz1,
                   mz2 = mz2, stringsAsFactors = FALSE)

  # iterate and calculate all combinatoins
  for(i in 1:nrow(df)) {

    # from adduct to neutral
    df$neutral1_from_adduct1[i] <- metabolomicsUtils::calc_neutral_mass(df$mz1[i], df$adduct1[i])
    df$neutral2_from_adduct2[i] <- metabolomicsUtils::calc_neutral_mass(df$mz2[i], df$adduct2[i])

    # from neutral to adduct
    df$adduct2_from_neutral1[i] <- metabolomicsUtils::calc_adduct_mass(df$neutral1_from_adduct1[i], df$adduct2[i])
    df$adduct1_from_neutral2[i] <- metabolomicsUtils::calc_adduct_mass(df$neutral2_from_adduct2[i], df$adduct1[i])

  }

  # select fitting adducts
  if(mzTolType == "abs") {

    filteredDf <- df[which(abs(df$mz1 - df$adduct1_from_neutral2) < mzTol & abs(df$mz2 - df$adduct2_from_neutral1) < mzTol),]

  } else if(mzTolType == "ppm") {

    filteredDf <- NULL

  } else {
    stop("unknown mzTolType")
  }

  if(!is.null(filteredDf) & nrow(filteredDf) > 0) {

    matchingResult <- ""

    for(i in 1:nrow(filteredDf)) {

      if(i == 1) {
        matchingResult <- paste0(filteredDf$adduct1[i], "<->", filteredDf$adduct2[i])
      } else {
        matchingResult <- paste0(matchingResult, " / ", filteredDf$adduct1[i], "<->", filteredDf$adduct2[i])
      }
    }

    # return resulting DF
    return(matchingResult)

  } else {

    return(NA)

  }
}

#' Checking masses from different ion modes
#'
#' This function is used to match two masses from the positive and negative ionization mode to check it they are potentially derived from the same metabolite. Based on the given adducts, all combinations are tested from positive to negative ion mode and the other way. If for both cases the mass error is below the given mass error it is assumed that the masses are derived from the same metabolite.
#'
#' @param posMz m/z value of positive ionization mode
#' @param negMz m/z value of negative ionization mode
#' @param posAdducts vector with adducts that are used for the positive ionization mode
#' @param negAdducts vector with adducts that are used for the negative ionization mode
#' @param mzTol m/z error, numeric value
#' @param mzTolType type of error used, absolute (abs) or relative (ppm)
#'
#' @examples
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{matchMassDiff}}
#' @export
match_mz_inter <- function(posMz, negMz, posAdducts, negAdducts, mzTol = 0.005, mzTolType = "abs") {

  # sanity checks
  # check if adduct definitions are good
  adduct_names <- metabolomicsUtils::get_adduct_names(mode = "all")

  if(!all(posAdducts %in% adduct_names)) {
    stop("one or more adducts in the posAdducts in not fitting")
  }

  if(!all(negAdducts %in% adduct_names)) {
    stop("one or more adducts in the negAdducts in not fitting")
  }

  # make combinations
  adductCombinations <- expand.grid(posAdducts, negAdducts)

  # make df
  df <- data.frame(posAdduct = as.character(adductCombinations$Var1),
                   negAdduct = as.character(adductCombinations$Var2),
                   posMz = posMz,
                   negMz = negMz, stringsAsFactors = FALSE)

  # iterate and calculate all combinatoins
  for(i in 1:nrow(df)) {

    # from pos to neg
    df$neutral_from_pos[i] <- metabolomicsUtils::calc_neutral_mass(df$posMz[i], df$posAdduct[i])
    df$neutral_from_neg[i] <- metabolomicsUtils::calc_neutral_mass(df$negMz[i], df$negAdduct[i])

    df$neg_from_pos[i] <- metabolomicsUtils::calc_adduct_mass(df$neutral_from_pos[i], df$negAdduct[i])
    df$pos_from_neg[i] <- metabolomicsUtils::calc_adduct_mass(df$neutral_from_neg[i], df$posAdduct[i])

  }

  # select fitting adducts
  if(mzTolType == "abs") {
    filteredDf <- df[which(abs(df$negMz - df$neg_from_pos) < mzTol & abs(df$posMz - df$pos_from_neg) < mzTol),]
  } else if(mzTolType == "ppm") {
    filteredDf <- NULL
  } else {
    stop("unknown mzTolType")
  }

  if(!is.null(filteredDf) & nrow(filteredDf) > 0) {
    matchingResult <- ""

    for(i in 1:nrow(filteredDf)) {

      if(i == 1) {
        matchingResult <- paste0(filteredDf$posAdduct[i], "<->", filteredDf$negAdduct[i])
      } else {
        matchingResult <- paste0(matchingResult, " / ", filteredDf$posAdduct[i], "<->", filteredDf$negAdduct[i])
      }
    }

    # return resulting DF
    return(matchingResult)
  } else {
    return(NA)
  }
}


#' Checking for mass diferences
#'
#' This functions is used to check of two measured masses are different by a give mass difference. This can be used to check if a neutral loss is found between two masses or if masses are connected by certain metabolic transformations, e.g. methylation, acetylation etc.
#'
#' @param mz1 First m/z value for which a mass difference shall be checked
#' @param mz2 Second m/z value for which a mass difference shall be checked
#' @param mzDiff Mass difference to check for
#' @param mzTol m/z error, numeric value
#' @param mzTolType type of error used, absolute (abs) or relative (ppm)
#'
#' @examples
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{matchMz}}
#' @export
matchMassDiff <- function(mz1, mz2, mzDiff, mzTol = 0.005, mzTolType = "abs") {

  match <- FALSE

  # check which error type
  if(mzTolType == "abs") {

    # check which mz is largerd
    if(mz1 > mz2) {

      if(abs(mz2 + mzDiff - mz1) < mzTol & abs(mz1 - mzDiff - mz2) < mzTol) {
        match <- TRUE
      }

    } else {

      if(abs(mz1 + mzDiff - mz2) < mzTol & abs(mz2 - mzDiff - mz1) < mzTol) {
        match <- TRUE
      }
    }

  } else if(mzTolType == "ppm") {

    # check which mz is largerd
    if(mz1 > mz2) {

      if(abs(mz2 + mzDiff - mz1) < (mzTol / 1e6 * mz1) & abs(mz1 - mzDiff - mz2) < (mzTol / 1e6 * mz1)) {
        match <- TRUE
      }

    } else {

      if(abs(mz1 + mzDiff - mz2) < (mzTol / 1e6 * mz2) & abs(mz2 - mzDiff - mz1) < (mzTol / 1e6 * mz2)) {
        match <- TRUE
      }
    }

  } else {
    stop("unknown mzTolType")
  }

  # return result
  return(match)
}

#' Calculate Kendrick Mass Defect
#'
#' This function calculates the Kendrick mass defect (KMD) for a given mass.
#'
#' @param mz Mass for which the KMD shall be calculated
#'
#' @example
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{calculateRefKendrickMassDefect}}
#' @seealso \code{\link{checkRmkd}}
#'
#' @export
calculateKendrickMassDefect <- function(mz) {

  # calculate KMD
  kendrickMass <- mz * 14 / 14.01565
  kmd <- kendrickMass %% 1

  # return Kendrick mass defect
  return(kmd)

}

#' Calculate the Kendrick Mass
#'
#' This function calculates the Kendrick mass for a give mass.
#'
#' @param mz Mass for which the KM shall be calculated
#'
#' @example
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{calculateRefKendrickMassDefect}}
#' @seealso \code{\link{calculateKendrickMassDefect}}
#' @seealso \code{\link{calculateNominalKendrickMass}}
#' @seealso \code{\link{checkRmkd}}
#'
#' @export
calculateKendrickMass <- function(mz) {

  # calculate KM
  kendrickMass <- mz * 14 / 14.01565

  # return Kendrick mass
  return(kendrickMass)
}

#' Calculate nominal Kendrick Mass
#'
#' This function calculates the nominal Kendrick Mass
#'
#' @param mz Mass for which the nominal KM shall be calculated
#'
#' @example
#' xxx
#'
#' @export
calculateNominalKendrickMass <- function(mz) {

  # calculate KM
  kendrickMass <- mz * 14 / 14.01565

  # nominal KM
  nominalKendrickMass <- as.integer(kendrickMass)

  # return nominal Kendrick Mass
  return(nominalKendrickMass)
}

#' Calculate referenced Kendrick Mass Defect
#'
#' This function calculates a referenced Kendrick Mass Defect (RKMD) for a given mass and a given reference KMD.
#'
#' @param mz Mass for which the RKMD shall be calculated
#' @param refkmd Reference Kendrick Mass Defect
#'
#' @examples
#' xxx
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{calculateKendrickMassDefect}}
#' @seealso \code{\link{checkRmkd}}
#'
#' @export
calculateRefKendrickMassDefect <- function(mz, refkmd) {

  # calculate KMD
  kmd <- calculateKendrickMassDefect(mz)

  # calculate RKMD
  rmkd <- (kmd - refkmd) / 0.013399

  # return referenced Kendrick mass defect
  return(rmkd)

}

#' Check if RMKD is negative integer
#'
#' This function checks if the RKMD is a negative integer
#'
#' @param rmkd calculated RMKD
#' @param error allowed error margin for RMKD
#'
#' @author Michael Witting, \email{michael.witting@@helmholtz-muenchen.de}
#'
#' @seealso \code{\link{calculateKendrickMassDefect}}
#' @seealso \code{\link{calculateRefKendrickMassDefect}}
#'
#' @export
checkRmkd <- function(rmkd, error = 0.15) {

  # rmkd can be only 0 or negative integer
  if(rmkd > 0) {
    return(FALSE)
  }

  # check if rmkd is in error range
  remainder <- rmkd %% -1

  if(abs(remainder) < error | abs(1 + remainder) < error) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
