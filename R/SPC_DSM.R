#' SPC_DSM
#' @description A single linear model with dispersion summation minimization
#' @details This function is to process phase error correction through a single linear model with dispersion summation minimization,
#' followed by Polynomial baseline correction
#' @param specdat A complex number vector of observed frequency domain data
#' @return A matrix with phase and baseline corrected spectra, the 1st column is for absorption spectrum, while the 2nd column is for dispersion spectrum.
#' @keywords 1D NMR, phase correction
#' @author Aixiang Jiang
#' @references
#' Chen, L., Weng, Z., Goh, L., & Garland, M. (2002). An efficient algorithm for automatic phase correction of NMR spectra based on
#' entropy minimization. Journal of Magnetic Resonance, 158, 1-2.
#'
#' Ernst, R. R. (1969). Numerical Hilbert transform and automatic phase correction in magnetic resonance spectroscopy.
#' Journal of Magnetic Resonance, 1, 7-26
#'
#' Kristian Hovde Liland, Trygve Almøy, Bjørn-Helge Mevik (2010), Optimal Choice of Baseline
#' Correction for Multivariate Calibration of Spectra, Applied Spectroscopy 64, pp. 1007-1016.
#' @import baseline
#' @export


SPC_DSM = function (specdat){

  ### the input is the real part of a spectrum
  ### output: phased real part of a spectrum

  hdat=cbind(Re(specdat), Im(specdat))

  ### set initial values as the same as my own methods
  #### find the max power spec value index
  pspec=hdat[,1]**2+hdat[,2]**2
  maxi=which.max(pspec) ### use the power spectrum is the same as using magnitude spectrum
  ph0Initial = -atan2(hdat[maxi,2],hdat[maxi,1]) ### it is in radian
  ph1Initial=0.005

  #### get optimized parameters of ph0 and ph1

  #### method1-4 has only diff in "fn" in the following line, all other lines are exactly the same
  optimRes=optim(par=c(ph0Initial,ph1Initial),fn=sumD, specDat=hdat)

  bestPh=optimRes$par

  #### now use this best ph to get phased data
  ####  -> first: get phase adding values
  nn=dim(hdat)[1]
  angles=bestPh[1]+bestPh[2]*c(1:nn)/nn

  ####  -> now, adding phase values
  dat3col=cbind(hdat, angles)
  phasedDat=t(apply(dat3col, 1, phaseCorr2)) ### output is a two column matrix: the phased real and the phased imaginary of freq data

  ##### return phased plus baseline corrected spectrum
  tryBL=baseline(t(phasedDat[,1]),method="modpolyfit")
  return(baseline::getCorrected(tryBL))

}
