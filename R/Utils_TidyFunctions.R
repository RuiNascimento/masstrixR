

# mutate wrapper fuction for rmkd
mutate_rmkd <- function(x, df) {

  # make grouping
  df %>% group_by(1:n())

  # iterate through all rmkds
  for(i in 1:length(x)) {

    rmkd <- x[i]
    varname1 <- paste0(names(rmkd))
    varname2 <- paste0(names(rmkd), "_check")

    df <- df %>%
      mutate()

  }
}

# mutate wrapper function to search for multiple neutral losses (x), pos mode
mutate_NL_pos <- function(x, df) {

  # make grouping
  df %>% group_by(1:n())

  # iterate through all neutral losses
  for(i in 1:length(x)) {

    nl <- x[i]
    varname <- paste0(names(nl), "_NL")

    df <- df %>%
      mutate(!!varname := matchMassDiff(Pos.mz, Pos.mz1, nl[[1]], mzTol = mzTol, mzTolType = mzTolType))
  }

  return(df)

}

# mutate wrapper function to search for multiple neutral losses (x), pos mode
mutate_NL_neg <- function(x, df) {

  # make grouping
  df %>% group_by(1:n())

  # iterate through all neutral losses
  for(i in 1:length(x)) {

    nl <- x[i]
    varname <- paste0(names(nl), "_NL")

    df <- df %>%
      mutate(!!varname := matchMassDiff(Neg.mz, Neg.mz1, nl[[1]], mzTol = mzTol, mzTolType = mzTolType))
  }

  return(df)

}
