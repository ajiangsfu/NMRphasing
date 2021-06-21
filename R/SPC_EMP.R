#' SPC_EMP
#' @description A single linear model with entropy minimization with negative peak penalty
#' @details This function is to process phase error correction through a single linear model with entropy minimization with negative peak penalty,
#' followed by Polynomial baseline correction
#' @param specdat A complex number vector of observed frequency domain data
#' @return A matrix with phase and baseline corrected spectra, the 1st column is for absorption spectrum, while the 2nd column is for dispersion spectrum.
#' @keywords 1D NMR, phase correction
#' @author Aixiang Jiang
#' @references
#' Binczyk F, Tarnawski R, Polanska J (2015) Strategies for optimizing the phase correction algorithms in Nuclear Magnetic Resonance spectroscopy. Biomed Eng Online 14 Suppl 2:S5.
#'
#' de Brouwer, H. (2009). Evaluation of algorithms for automated phase correction of NMR spectra. J Magn Reson, 201, 230-238.
#'
#' Liland KH, Almøy T, Mevik B (2010), Optimal Choice of Baseline
#' Correction for Multivariate Calibration of Spectra, Applied Spectroscopy 64, pp. 1007-1016.
#' @import baseline
#' @export


SPC_EMP = function (specdat){
  hdat=cbind(Re(specdat), Im(specdat))

  pspec=hdat[,1]**2+hdat[,2]**2
  maxi=which.max(pspec)
  ph0Initial = -atan2(hdat[maxi,2],hdat[maxi,1])
  ph1Initial=0.005

  #### get optimized parameters of ph0 and ph1

  optimRes=optim(par=c(ph0Initial,ph1Initial),fn=entropyP, specDat=hdat)
  bestPh=optimRes$par

  nn=dim(hdat)[1]
  angles=bestPh[1]+bestPh[2]*c(1:nn)/nn

  dat3col=cbind(hdat, angles)
  phasedDat=t(apply(dat3col, 1, phaseCorr2)) ### output is a two column matrix: the phased real and the phased imaginary of freq data

  ##### return phased plus baseline corrected spectrum
  tryBL=baseline(t(phasedDat[,1]),method="modpolyfit")
  return(baseline::getCorrected(tryBL)[1,])
}
