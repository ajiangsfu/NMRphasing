#' A function rather aimed at developers
#' @import MassSpecWavelet
#' @noRd

peakSearch=function(datin){ ### datin is a vector
  peakInfo = MassSpecWavelet::peakDetectionCWT(datin)
  majorPeakInfo = peakInfo$majorPeakInfo
  peakIndex = majorPeakInfo$peakIndex
  return(peakIndex)
}
