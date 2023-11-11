#' SPC_DSM
#' @description A single linear model with dispersion summation minimization.
#' @details This function is to process phase error correction through a single linear model with dispersion summation minimization,
#' followed by Polynomial baseline correction.
#' @param specdat A complex number vector of observed frequency domain data
#' @return A numeric vector of phase corrected absorption spectrum
#' @concept phase correction
#' @author Aixiang Jiang
#' @references
#' Chen, L., Weng, Z., Goh, L., & Garland, M. (2002). An efficient algorithm for automatic phase correction of NMR spectra based on
#' entropy minimization. Journal of Magnetic Resonance, 158, 1-2.
#'
#' Ernst, R. R. (1969). Numerical Hilbert transform and automatic phase correction in magnetic resonance spectroscopy.
#' Journal of Magnetic Resonance, 1, 7-26
#'
#' Liland KH, Almøy T, Mevik B (2010), Optimal Choice of Baseline
#' Correction for Multivariate Calibration of Spectra, Applied Spectroscopy 64, pp. 1007-1016.

#' @import baseline
#'
#' @examples
#' data("fdat")
#' spc_dsm_phased1 <- SPC_DSM(fdat$frequency_domain)

#' @export


SPC_DSM = function (specdat){

  hdat=cbind(Re(specdat), Im(specdat))

  pspec=hdat[,1]**2+hdat[,2]**2
  maxi=which.max(pspec)
  ph0Initial = -atan2(hdat[maxi,2],hdat[maxi,1])
  ph1Initial=0.005

  #### get optimized parameters of ph0 and ph1
  optimRes=stats::optim(par=c(ph0Initial,ph1Initial),fn=sumD, specDat=hdat)

  bestPh=optimRes$par

  nn=dim(hdat)[1]
  angles=bestPh[1]+bestPh[2]*c(1:nn)/nn

  dat3col=cbind(hdat, angles)
  phasedDat=t(apply(dat3col, 1, phaseCorr2))

  ##### return phased plus baseline corrected spectrum
  tryBL=baseline(t(phasedDat[,1]),method="modpolyfit")
  return(baseline::getCorrected(tryBL)[1,])

}





