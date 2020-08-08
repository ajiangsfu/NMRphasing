#' A function rather aimed at developers
#' @noRd

EMP = function (specIn){
  hdat=cbind(Re(specIn), Im(specIn))
  ### set initial values as the same as my own methods
  #### find the max power spec value index
  pspec=hdat[,1]**2+hdat[,2]**2
  maxi=which.max(pspec) ### use the power spectrum is the same as using magnitude spectrum
  ph0Initial = -atan2(hdat[maxi,2],hdat[maxi,1]) ### it is in radian
  ph1Initial=0.005

  #### get optimized parameters of ph0 and ph1

  optimRes=optim(par=c(ph0Initial,ph1Initial),fn=entropyP, specIn=hdat)
  bestPh=optimRes$par

  #### now use this best ph to get phased data
  ####  -> first: get phase adding values
  nn=dim(hdat)[1]
  angles=bestPh[1]+bestPh[2]*c(1:nn)/nn

  ####  -> now, adding phase values
  dat3col=cbind(hdat, angles)
  phasedDat=t(apply(dat3col, 1, phaseCorr2)) ### output is a two column matrix: the phased real and the phased imaginary of freq data

  return(phasedDat[,1])
}
