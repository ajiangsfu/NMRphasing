#' A function primarily aimed at developers
#' The original 'localMaximumSlidingWindow' from the MassSpecWavelet package has been revised and converted to 'localMaximumSlidingWindowR'.
#' The conversion was made to address specific bugs and ensure compatibility with the current workflow.
#' The function findLocalMaxWinSizeR is adapted from 'find_local_maximum.c' in the MassSpecWavelet package
#' @noRd

localMaximumR <- function(x, winSize = 5) {
  algo <- getOption("MassSpecWavelet.localMaximum.algorithm", "faster")
  if (!algo %in% c("new", "classic", "faster")) {
    warning('Invalid algorithm "', algo, '". Use either "new", "faster" or "classic". Assuming "faster".')
    algo <- "faster"
  }
  if (algo == "new") {
    local_max <- findLocalMaxWinSizeR(x, capWinSize = winSize)
      # since findLocalMaxWinSize is not included in MassSpecWavelet's NAMESPACE, use :::
    localMax <- as.integer(local_max >= winSize)
    return(localMax)
  } else if (algo == "faster") {
    localMax <- localMaximumSlidingWindowR(x, winSize)
    ## Check whether there is some local maxima have in between distance less than winSize
    maxInd <- which(localMax > 0)
    selInd <- which(diff(maxInd) < winSize)
    if (length(selInd) > 0) {
      selMaxInd1 <- maxInd[selInd]
      selMaxInd2 <- maxInd[selInd + 1L]
      temp <- x[selMaxInd1] - x[selMaxInd2]
      localMax[selMaxInd1[temp <= 0]] <- 0L
      localMax[selMaxInd2[temp > 0]] <- 0L
    }
    return(localMax)
  }
  len <- length(x)
  rNum <- ceiling(len / winSize)

  ## Transform the vector as a matrix with column length equals winSize
  ## and find the maximum position at each row.
  y <- matrix(c(x, rep(x[len], rNum * winSize - len)), nrow = winSize)
  y.maxInd <- apply(y, 2, which.max)
  ## Only keep the maximum value larger than the boundary values
  selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))

  ## keep the result
  localMax <- rep(0L, len)
  localMax[(selInd - 1) * winSize + y.maxInd[selInd]] <- 1L

  ## Shift the vector with winSize/2 and do the same operation
  shift <- floor(winSize / 2)
  rNum <- ceiling((len + shift) / winSize)
  y <- matrix(c(rep(x[1], shift), x, rep(x[len], rNum * winSize - len - shift)), nrow = winSize)
  y.maxInd <- apply(y, 2, which.max)
  ## Only keep the maximum value larger than the boundary values
  selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
  localMax[(selInd - 1) * winSize + y.maxInd[selInd] - shift] <- 1L

  ## Check whether there is some local maxima have in between distance less than winSize
  maxInd <- which(localMax > 0)
  selInd <- which(diff(maxInd) < winSize)
  if (length(selInd) > 0) {
    selMaxInd1 <- maxInd[selInd]
    selMaxInd2 <- maxInd[selInd + 1L]
    temp <- x[selMaxInd1] - x[selMaxInd2]
    localMax[selMaxInd1[temp <= 0]] <- 0L
    localMax[selMaxInd2[temp > 0]] <- 0L
  }

  return(localMax)
}
