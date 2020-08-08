#' MPC_EMP
#' @description Multiple single linear models with entropy minimization with negative peak penalty
#' @details This function is to process phase error correction through multiple single linear models with entropy minimization with negative peak penalty,
#' followed by Polynomial baseline correction
#' @param specdat A complex number vector of observed frequency domain data
#' @return A matrix with phase and baseline corrected spectra, the 1st column is for absorption spectrum, while the 2nd column is for dispersion spectrum.
#' @keywords 1D NMR, phase correction
#' @author Aixiang Jiang
#' @references
#'
#' This is our new algorithm
#'
#' de Brouwer, H. (2009). Evaluation of algorithms for automated phase correction of NMR spectra. J Magn Reson, 201, 230-238.
#'
#' Kristian Hovde Liland, Trygve Almøy, Bjørn-Helge Mevik (2010), Optimal Choice of Baseline
#' Correction for Multivariate Calibration of Spectra, Applied Spectroscopy 64, pp. 1007-1016.
#' @export


MPC_EMP = function(specdat){

  cplxDat=specdat
  pp=(Re(cplxDat))**2+(Im(cplxDat))**2
  mm=sqrt(pp)

  peakIndex=peakSearch(mm)
  k=length(peakIndex)

  peak1=peakIndex[-k]
  peak2=peakIndex[-1]
  peaks=cbind(peak1,peak2)

  #### find valley index
  valleyIndex=apply(peaks,1,FUN=function(x){
    mdat=mm[x[1]:x[2]]
    which.min(mdat)+x[1]-1
  })

  nn=length(specdat)
  valleyIndex=c(1,valleyIndex, nn)

  vv=length(valleyIndex)
  valleyL=valleyIndex[1:(vv-1)]
  valleyL[2:(vv-1)]=valleyL[2:(vv-1)]+1
  valleyR=valleyIndex[2:vv]
  valleys=cbind(valleyL,valleyR)

  phasedComb=apply(valleys,1,FUN=function(x){
    res=EMP(cplxDat[x[1]:x[2]])
  })

  phasedAll=unlist(phasedComb)
  tryBL=myBaseline(phasedAll,bsDf=5, BL_method="modpolyfit")
  return(tryBL)

}
