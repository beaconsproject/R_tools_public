ks_plot <- function(refVal, netVal, plotTitle="") {
  
  regLab <- "Reference area"
  netLab <- "Network"
  
  z1 <- c(refVal, netVal)
  z2 <- c(rep(regLab,length(refVal)), rep(netLab,length(netVal)))
  zz <- data.frame(cbind(z1,z2),stringsAsFactors=FALSE)
  names(zz) <- c("values","criteria")
  zz$values <- round(as.numeric(zz$values),3)
  
  # create and save density plot
  p <- ggplot2::ggplot(zz, ggplot2::aes(x=.data$values)) + ggplot2::geom_density(ggplot2::aes(group=.data$criteria, color=.data$criteria)) +
    ggplot2::ggtitle(plotTitle) +
    ggplot2::labs(x="Indicator value", y="Density")
  
  return(p)
}
