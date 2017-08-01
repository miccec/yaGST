plot.mwwGST <-
function (obj, rankedList = NULL, main = "MWW-GST", ...) {
  if (length(obj) < 2) {
    message("Nothing to plot.")
    return(invisible(NULL))
  }
  if(is.null(rankedList) & is.null(obj$rankedList)) {
    message("The ranked-list is needed!")
    message("Nothing to plot.")
    return(invisible(NULL))
  }
  
  if(is.null(rankedList)) rankedList <- obj$rankedList
    
  require(ggplot2)
  
  K <- 500
  cValue <- quantile(rankedList, probs = 0:K/K)
  tmp <- unique(cValue)
  if(length(cValue) != length(tmp)) {
    message("The nr of knots is less than 500.")
    cValue <- tmp
  }
  m_J <- obj$actualGeneSetCount
  F_J <- rankedList[obj$actualGeneSet]
  F_J <- cumsum(table(cut(F_J, breaks = cValue, include.lowest = TRUE)))/m_J
  m_Jc <- obj$lengthOfRankedList - m_J
  F_Jc <- rankedList[setdiff(names(rankedList), obj$actualGeneSet)]
  F_Jc <- cumsum(table(cut(F_Jc, breaks = cValue, include.lowest = TRUE)))/m_Jc
  ddata <- data.frame(F_Jc = F_Jc, F_J = F_J)
  
  ggp <- ggplot(ddata, aes(x = F_Jc, y= F_J)) 
  ggp <- ggp + labs(title=main, y = "gene set distribution function",
                    x = "genes distribution function")
  ggp <- ggp + geom_line(col = "red", size = 2)
  ggp <- ggp + geom_abline(intercept = 0, slope = 1)
  tmp <- paste0("NES = ", round(obj$nes, 4), "; pValue = ", 
                format(obj$p.value, digit = 4, scientific = TRUE))
  ggp <- ggp + annotate("text", x = 0.02, y = 0.98, label = tmp, size = 7, adj = 0)
  if(!is.null(obj$geneSetName)) 
    ggp <- ggp + annotate("text", x = 0.02, y = 1.02, label = obj$geneSetName, size = 7, adj = 0)
  ggp <- ggp + geom_ribbon(aes(ymin=F_J, ymax=1, x=F_Jc, fill = "band"), alpha = 0.2)
  ggp <- ggp + theme(legend.position = "none")
  return(ggp)
}
