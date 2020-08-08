#' NMRphasing
#' @description Phase error correction wrap up function
#' @details This is a wrap function to process phase error correction with seven different choices, followed by Polynomial baseline correction when necessary
#' @param specDatIn Input spectrum data, which can be one of the four formats:
#'                a vector of absorption specrtrum; a complex vector; a data matrix or a data frame with two columns of spectrum data,
#'                which 1st column is for absorption spectrum, and 2nd column is for dispersion spectrum
#' @param absorptionOnly A logical variable to tell us if specDatIn is a a vector of absorption specrtrum, default is false
#' @param method One of phase correction method. There are seven methods right now, which are "NLS", "MPC_DAOM", "MPC_EMP", "SPC_DAOM", "SPC_EMP", "SPC_AAM", "SPC_DSM",
#'               with NLS, non-linear shrinkage as default
#' @return A matrix with phase and baseline corrected spectra, the 1st column is for absorption spectrum, while the 2nd column is for dispersion spectrum.
#' @keywords 1D NMR, phase correction, a single linear model, minimization of absolute area
#' @author Aixiang Jiang
#' @references
#' de Brouwer, H. (2009). Evaluation of algorithms for automated phase correction of NMR spectra. J Magn Reson, 201, 230-238.
#'
#' Dzakula, Z. (2000). Phase angle measurement from peak areas (PAMPAS). J Magn Reson, 146, 20-32.
#'
#' Chen, L., Weng, Z., Goh, L., & Garland, M. (2002). An efficient algorithm for automatic phase correction of NMR spectra based on
#' entropy minimization. Journal of Magnetic Resonance, 158, 1-2.
#'
#' Ernst, R. R. (1969). Numerical Hilbert transform and automatic phase correction in magnetic resonance spectroscopy.
#' Journal of Magnetic Resonance, 1, 7-26
#'
#' Kristian Hovde Liland, Trygve Almøy, Bjørn-Helge Mevik (2010), Optimal Choice of Baseline
#' Correction for Multivariate Calibration of Spectra, Applied Spectroscopy 64, pp. 1007-1016.
#'
#'
#' @export


NMRphasing = function (specDatIn, absorptionOnly = FALSE, method = c("NLS", "MPC_DAOM", "MPC_EMP", "SPC_DAOM", "SPC_EMP", "SPC_AAM", "SPC_DSM")){
  ## this function can accept three formats of specDatIn
  ## specDatIn can be a vector of absorption specrtrum, a complex vector, or a data matrix/frame with two columns of spectrum data
  ##       if it is a data frame or data matrix, the 1st column for absorption, and the 2nd column is for dispersion
  ## right now, there are only 7 fixed methods, in the future, should allow users to choose:
  ##  i) single linear model, or multiple linear models, or NLS (shrinkage) that can not be combined with the 2) step though
  ##  ii) one of optimization functions, right now, only 4 of them, should add more
  ##  iii) detect baseline bias and then use polynomial methods to correct it, or without detection step
  ##           or maybe add more choices including no baseline correction choices in the future
  datin = NA
  if(absorptionOnly){
    datin = HilbertWithFT(specDatIn)
  }else if(is.complex(specDatIn)){
    datin = specDatIn
  }else{
    datin = complex(real = specDatIn[,1], imaginary = specDatIn[,2])
  }

  outdat = NA

  ### the following part is for choices of 7 methods, all of them should have a output object: outdat
  if(method == "MPC_DAOM"){
    outdat = MPC_DAOM(specdat = datin)
  }else if(method == "MPC_EMP"){
     outdat = MPC_EMP(specdat = datin)
  }else if(method == "SPC_DAOM"){
     outdat = SPC_DAOM(specdat = datin)
  }else if(method == "SPC_EMP"){
     outdat = SPC_EMP(specdat = datin)
  }else if(method == "SPC_AAM"){
     outdat = SPC_AAM(specdat = datin)
  }else if(method == "SPC_DSM"){
     outdat = SPC_DSM(specdat = datin)
  }else{  ## default is NLS
     outdat = NLS(specdat = datin)
  }

  return(outdat)

}
