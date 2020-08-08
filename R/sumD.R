#' A function rather aimed at developers
#' @noRd

sumD = function(phasePara, specDat) {

  #### phasePara contains a pair of phase adding parameter for a linear model, they are in radius
  ####   phasePara[1] is for intercept, and phasePara[2] is for slope
  #### specDat is a matrix with two columns: the real and the imaginary parts of frequency domain data

  #### when set phasePara in optim function, set as: optim(par=c(ph0Initial,ph1Initial),fn=absArea, specDat=...)
  ####   of course, I need to define ph0Initial and ph1Initial before calling optim function

  #### use a linear model to get phase adding values for a whole spectrum based on phasePara
  #### get length of a spectrum
  n=dim(specDat)[1]
  phases=phasePara[1]+phasePara[2]*c(1:n)/n

  dat3col=cbind(specDat, phases)
  phasedDat=t(apply(dat3col, 1, phaseCorr2)) ### output is a two column matrix: the phased real and the phased imaginary of freq data
  return(sum(phasedDat[,2]))
}
