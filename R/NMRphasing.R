#' NMRphasing
#' @description Phase error correction wrap up function
#' @details This is a wrap function to process phase error correction and baseline correction with nine different choices.
#' @param specDatIn Input spectrum data, which can be one of the four formats:
#'                a vector of absorption spectrum; a complex vector; a data matrix or a data frame with two columns of spectrum data,
#'                which 1st column is for absorption spectrum, and 2nd column is for dispersion spectrum
#' @param absorptionOnly A logical variable to tell us if specDatIn is a a vector of absorption specrtrum, default is false
#' @param method One of phase correction and baseline correction methods. There are nine available choices, which are "NLS", "MPC_DAOM", "MPC_EMP", "MPC_AAM", "MPC_DSM", "SPC_DAOM", "SPC_EMP", "SPC_AAM", "SPC_DSM",
#'               with "NLS", non-linear shrinkage as default.
#' @return A numeric vector of phase corrected absorption spectrum
#' @concept phase error correction
#' @author Aixiang Jiang
#' @references

#' Binczyk F, Tarnawski R, Polanska J (2015) Strategies for optimizing the phase correction algorithms in Nuclear Magnetic Resonance spectroscopy. Biomed Eng Online 14 Suppl 2:S5.
#'
#' Chen L, Weng Z, Goh L, Garland M (2002) An efficient algorithm for automatic phase correction of NMR spectra based on entropy minimization. J Magn Reson 158:164–168.
#'
#' de Brouwer H (2009) Evaluation of algorithms for automated phase correction of NMR spectra. J Magn Reson 201:230–238.
#'
#' Džakula Ž (2000) Phase Angle Measurement from Peak Areas (PAMPAS). J Magn Reson 146:20–32.
#'
#' Ernst RR (1969) Numerical Hilbert transform and automatic phase correction in magnetic resonance spectroscopy. J Magn Reson 1969 1:7–26.
#'
#' Liland KH, Almøy T, Mevik B (2010), Optimal Choice of Baseline
#' Correction for Multivariate Calibration of Spectra, Applied Spectroscopy 64, pp. 1007-1016.

#' @examples
#' data("fdat")
#' nls_phased <- NMRphasing(specDatIn = fdat$frequency_domain, method = "NLS")


#' @export


NMRphasing = function (specDatIn, absorptionOnly = FALSE,
                       method = c("NLS", "MPC_DANM", "MPC_EMP","MPC_AAM", "MPC_DSM","SPC_DANM", "SPC_EMP", "SPC_AAM", "SPC_DSM")){
  datin = NA
  if(absorptionOnly){
    datin = HilbertWithFT(specDatIn)
  }else if(is.complex(specDatIn)){
    datin = specDatIn
  }else{
    datin = complex(real = specDatIn[,1], imaginary = specDatIn[,2])
  }

  outdat = NA

  if(method == "MPC_DANM"){
     outdat = MPC_DANM(specdat = datin)
  }else if(method == "MPC_EMP"){
     outdat = MPC_EMP(specdat = datin)
  }else if(method == "MPC_AAM"){
    outdat = MPC_AAM(specdat = datin)
  }else if(method == "MPC_DSM"){
    outdat = MPC_DSM(specdat = datin)
  }else if(method == "SPC_DANM"){
     outdat = SPC_DANM(specdat = datin)
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

