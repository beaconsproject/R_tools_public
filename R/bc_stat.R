bc_stat <- function(refVal, netVal) {
  
  x1 <- dplyr::as_tibble(refVal) %>%
    dplyr::count(.data$value)
  names(x1) <- c("cat","strata")
  
  x2 <- dplyr::as_tibble(netVal) %>%
    dplyr::count(.data$value)
  names(x2) <- c("cat","reserve")
  
  x <- merge(x1,x2,by="cat",all=T)
  #x$cat <- as.character(x$cat)
  x$strata <- as.numeric(x$strata)
  x$reserve <- as.numeric(x$reserve)
  x$reserve[is.na(x$reserve)] <- 0
  x$reserve[is.na(x$strata)] <- 0 # this is needed in case there is one reserve pixel and no strata pixel
  x$strata[is.na(x$strata)] <- 0
  x$strata <- x$strata/sum(x$strata)
  x$reserve <- x$reserve/sum(x$reserve)
  
  # calculate Bray-Curtis dissimilariy
  ri <- round(sum(abs(x$strata-x$reserve))/(sum(x$strata)+sum(x$reserve)), 3)
  return(ri)
}
