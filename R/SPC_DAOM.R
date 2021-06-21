
#' SPC_DAOM
#' @description A single linear model with Minimization of difference between absolute area and ordinary area
#' @details This function is to process phase error correction through a single linear model with Minimization of difference between absolute area and ordinary area,
#' followed by Polynomial baseline correction when necessary
#' @param specdat A complex number vector of observed frequency domain data
#' @return A matrix with phase and baseline corrected spectra, the 1st column is for absorption spectrum, while the 2nd column is for dispersion spectrum.
#' @keywords 1D NMR, phase correction
#' @author Aixiang Jiang
#' @references
#' Jiang A, Hanley JA, and Nadon R, 2021, 1D NMR Phase Error Correction with New Modeling Methods (in preparation)
#'
#' Jiang A, Gravel A, Hanley JA, and Nadon R, 2021, Comparing phase error correction methods with NMR spike-in experiments (in preparation)
#'
#' Liland KH, Almøy T, Mevik B (2010), Optimal Choice of Baseline
#' Correction for Multivariate Calibration of Spectra, Applied Spectroscopy 64, pp. 1007-1016.
#' @export


SPC_DAOM =function (specdat){

  hdat=cbind(Re(specdat), Im(specdat))

  pspec=hdat[,1]**2+hdat[,2]**2
  maxi=which.max(pspec)
  ph0Initial = -atan2(hdat[maxi,2],hdat[maxi,1])
  ph1Initial=0.005

  #### get optimized parameters of ph0 and ph1
  optimRes=optim(par=c(ph0Initial,ph1Initial),fn=areaDiff, specDat=hdat)
  bestPh=optimRes$par

  nn=dim(hdat)[1]
  angles=bestPh[1]+bestPh[2]*c(1:nn)/nn

  dat3col=cbind(hdat, angles)
  phasedDat=t(apply(dat3col, 1, phaseCorr2)) ### output is a two column matrix: the phased real and the phased imaginary of freq data

  tryBL=myBaseline(phasedDat[,1],bsDf=5, BL_method="modpolyfit")
  return(tryBL)

}
