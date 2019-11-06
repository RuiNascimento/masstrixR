# plotting function
#'
#'
#' @import ggplot2
#' @import grid
#' @export
makeMirrorPlot <- function(x, y, align = FALSE, plotIt = FALSE, savePlot = FALSE, fileName = "plot.png", mzTol = 0.005, treshold = 0.01, title = "Mirrorplot", xlim = NULL, ...) {


  # no further continued
  .Deprecated(msg = "'searchByPrecurosr' will be removed in the next version")


  if(align) {

    #use aligned spectra
    alignedSpectra <- alignSpectra(x, y, mzTol = mzTol)
    commonPeaks <- alignedSpectra[which(alignedSpectra$intensity.top > 0 & alignedSpectra$intensity.bottom > 0),]

    noPeaks_x <- length(mz(x))
    noPeaks_y <- length(mz(y))

    # check xlim
    if(is.null(xlim)) {
      xlim <- c(min(alignedSpectra$mz) - 5, max(alignedSpectra$mz) + 5)
    }

    if(nrow(commonPeaks) > 0) {

      p1 <- ggplot() +

        # upper spectrum
        geom_linerange(data = alignedSpectra, aes(x = mz, ymin = 0, ymax = intensity.top), colour = "blue") +
        annotate("text", x = min(alignedSpectra$mz) - 5, y = 120, label = paste0("Precursor m/z: ", round(precursorMz(x), 4), " / Common Peaks: ", nrow(commonPeaks), " / Total Peaks: ", noPeaks_x), hjust = 0) +

        # lower spectrum
        geom_linerange(data = alignedSpectra, aes(x = mz, ymax = 0, ymin = intensity.bottom * -1), colour = "red") +
        annotate("text", x = min(alignedSpectra$mz) - 5, y = -120, label = paste0("Precursor m/z: ", round(precursorMz(y), 4), " / Common Peaks: ", nrow(commonPeaks), " / Total Peaks: ", noPeaks_y), hjust = 0) +

        # common peaks marker
        geom_point(data = commonPeaks, aes(x = mz, y = intensity.top + 5), shape = 25, colour = "black", fill = "blue") +
        geom_point(data = commonPeaks, aes(x = mz, y = intensity.bottom * - 1 -5), shape = 24, colour = "black", fill = "red") +

        # title
        ggtitle(title) +

        # scaling
        scale_x_continuous(limits = xlim) +
        scale_y_continuous(breaks = c(-100, -50, 0, 50, 100)) +

        # axis
        xlab("m/z") + ylab("normalized intensity") +

        # theme
        theme_bw() +
        theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black", linetype = "solid"))

      if(plotIt) {
        plot(p1)
      }

      if(savePlot) {
        print(paste0("saving plot to: ", fileName))
        ggsave(fileName, plot = p1, width = 150, height = 75, units = "mm", dpi = 300, limitsize = FALSE)
      }


    } else {

      p1 <- ggplot() +

        # upper spectrum
        geom_linerange(data = alignedSpectra, aes(x = mz, ymin = 0, ymax = intensity.top), colour = "blue") +
        annotate("text", x = min(alignedSpectra$mz) - 5, y = 120, label = paste0("Precursor m/z: ", round(precursorMz(x), 4), " / Common Peaks: ", nrow(commonPeaks)), hjust = 0) +

        # lower spectrum
        geom_linerange(data = alignedSpectra, aes(x = mz, ymax = 0, ymin = intensity.bottom * -1), colour = "red") +
        annotate("text", x = min(alignedSpectra$mz) - 5, y = -120, label = paste0("Precursor m/z: ", round(precursorMz(y), 4), " / Common Peaks: ", nrow(commonPeaks)), hjust = 0) +

        # title
        ggtitle(title) +

        # scaling
        scale_x_continuous(limits = xlim) +
        scale_y_continuous(breaks = c(-100, -50, 0, 50, 100)) +

        # axis
        xlab("m/z") + ylab("normalized intensity") +

        # theme
        theme_bw() +
        theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black", linetype = "solid"))

      if(plotIt) {
        plot(p1)
      }

      if(savePlot) {
        print(paste0("saving plot to: ", fileName))
        ggsave(fileName, plot = p1, width = 150, height = 75, units = "mm", dpi = 300, limitsize = FALSE)
      }

    }
  } else {
    NULL
  }
}

# plotting function
#'
#'
#'
#' @import ggplot2
#' @import grid
#' @export
plotSpectrum <- function(x, highlight = FALSE, highlightMz = NULL, mzTol = 0.005, plotIt = FALSE, xlim = NULL, ...) {

  # check for class
  if(class(x) == "Spectra") {
    x <- x[[1]]
  }


  if(highlight) {

    #check if m/z values for highlighting are supplied
    if(is.null(highlightMz)) {
      stop("No m/z to highlight supplied!")
    }

    # make spectrum object for highlighting
    highlightSpectrum <- new("Spectrum1",
                             mz = highlightMz,
                             intensity = rep(1, length(highlightMz)))

    #use aligned spectra
    alignedSpectra <- alignSpectra(highlightSpectrum, x, mzTol = mzTol)
    commonPeaks <- alignedSpectra[which(alignedSpectra$intensity.top > 0 & alignedSpectra$intensity.bottom > 0),]

    # check xlim
    if(is.null(xlim)) {
      xlim <- c(min(alignedSpectra$mz) - 5, max(alignedSpectra$mz) + 5)
    }

    if(nrow(commonPeaks) > 0) {

      p1 <- ggplot() +

        # upper spectrum
        geom_linerange(data = alignedSpectra, aes(x = mz, ymin = 0, ymax = intensity.bottom), colour = "blue") +

        # common peaks marker
        geom_point(data = commonPeaks, aes(x = mz, y = intensity.bottom + 5), shape = 25, colour = "black", fill = "blue") +

        # title
        ggtitle(paste0("Precursor m/z: ", round(precursorMz(x), 4), " / Common Peaks: ", nrow(commonPeaks))) +

        # scaling
        scale_x_continuous(limits = xlim) +
        scale_y_continuous(breaks = c(0, 50, 100)) +

        # axis
        xlab("m/z") + ylab("normalized intensity") +

        # theme
        theme_bw() +
        theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black", linetype = "solid"))

      if(plotIt) {
        plot(p1)
      }

    } else {

      p1 <- ggplot() +

        # upper spectrum
        geom_linerange(data = alignedSpectra, aes(x = mz, ymin = 0, ymax = intensity.top), colour = "blue") +

        # title
        ggtitle(paste0("Precursor m/z: ", round(precursorMz(x), 4))) +

        # scaling
        scale_x_continuous(limits = xlim) +
        scale_y_continuous(breaks = c(0, 50, 100)) +

        # axis
        xlab("m/z") + ylab("normalized intensity") +

        # theme
        theme_bw() +
        theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black", linetype = "solid"))

      if(plotIt) {
        plot(p1)
      }

    }
  } else {

    alignedSpectra <- data.frame(mz = mz(x),
                                 intensity.top = intensity(x) / max(intensity(x) * 100))

    p1 <- ggplot() +

      # upper spectrum
      geom_linerange(data = alignedSpectra, aes(x = mz, ymin = 0, ymax = intensity.top), colour = "blue") +

      # title
      ggtitle(paste0("Precursor m/z: ", round(precursorMz(x), 4))) +

      # scaling
      scale_x_continuous(limits = xlim) +
      scale_y_continuous(breaks = c(0, 50, 100)) +

      # axis
      xlab("m/z") + ylab("normalized intensity") +

      # theme
      theme_bw() +
      theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black", linetype = "solid"))

    if(plotIt) {
      plot(p1)
    }

  }
}


# plotting function
#'
#'
#' @import ggplot2
#' @import grid
#' @export
make_mirror_plot <- function(x, y, align = FALSE, savePlot = FALSE,
                             fileName = "plot.png", returnPlot = FALSE,
                             mzTol = 0.005, treshold = 0.01,
                             title = "Mirrorplot", xlim = NULL, ...) {

  # check if spectra shall be aligned or binned
  if(align) {

    #use aligned spectra
    compared_spectra <- alignSpectra(x, y, mzTol = mzTol, mzTolType = mzTolType)

  } else {

    # use binned spectra
    compared_spectra <- bin_Spectra(x, y, ...)

  }

  # get number of total peaks in each spectrum
  noPeaks_x <- length(mz(x))
  noPeaks_y <- length(mz(y))

  # check xlim
  if(is.null(xlim)) {
    xlim <- c(min(compared_spectra$mz) - 5, max(compared_spectra$mz) + 5)
  }

  # create annotation for spectra
  upper_annotation <- paste0("Precursor m/z: ", round(precursorMz(x), 4), " / Common Peaks: ", nrow(commonPeaks), " / Total Peaks: ", noPeaks_x)
  lower_annotation <- paste0("Precursor m/z: ", round(precursorMz(y), 4), " / Common Peaks: ", nrow(commonPeaks), " / Total Peaks: ", noPeaks_y)

  # create base plot
  p1 <- ggplot() +

    # upper spectrum
    geom_linerange(data = compared_spectra, aes(x = mz, ymin = 0, ymax = intensity.top), colour = "blue") +
    annotate("text", x = min(compared_spectra$mz) - 5, y = 120, label = upper_annotation, hjust = 0) +

    # lower spectrum
    geom_linerange(data = compared_spectra, aes(x = mz, ymax = 0, ymin = intensity.bottom * -1), colour = "red") +
    annotate("text", x = min(compared_spectra$mz) - 5, y = -120, label = , hjust = 0) +

    # title
    ggtitle(title) +

    # scaling
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(breaks = c(-100, -50, 0, 50, 100)) +

    # axis
    xlab("m/z") + ylab("normalized intensity") +

    # theme
    theme_bw() +
    theme(panel.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", linetype = "solid"))

  # if common peak exist highlight them
  if(nrow(commonPeaks) > 0) {
    p1 <- p1 +     # common peaks marker
      geom_point(data = commonPeaks, aes(x = mz, y = intensity.top + 5), shape = 25, colour = "black", fill = "blue") +
      geom_point(data = commonPeaks, aes(x = mz, y = intensity.bottom * - 1 -5), shape = 24, colour = "black", fill = "red")
  }

  # plot the plot
  plot(p1)

  # plot saving options
  if(savePlot) {
    print(paste0("saving plot to: ", fileName))
    ggsave(fileName, plot = p1, width = 150, height = 75, units = "mm", dpi = 300, limitsize = FALSE)
  }

  # check if plot shall be returned
  if(returnPlot) {
    return(p1)
  }
}
