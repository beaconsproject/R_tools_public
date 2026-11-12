# This should be written in RCPP
ks_stat <- function(refVal, netVal) {
  # calculate KS statistic (representation index)
  ri <- suppressWarnings(round(ks.test(refVal, netVal)[[1]][[1]], 3))
  return(ri)
}
