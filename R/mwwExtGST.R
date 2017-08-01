mwwExtGST <-
function(rankedList, geneSetUp, geneSetDown, minLenGeneSet = 15, moreDetails = FALSE, verbose = TRUE) {
  ccall <- match.call()
  doubleRankedList <- c(rankedList, -rankedList)
  flag <- c(rep(TRUE, length(rankedList)), rep(FALSE, length(rankedList)))
  oorder <- order(doubleRankedList, decreasing = TRUE)
  doubleRankedList <- doubleRankedList[oorder]
  flag <- flag[oorder]
  hits <- (names(doubleRankedList) %in% geneSetUp) & flag
  hits <- hits | ((names(doubleRankedList) %in% geneSetDown) & !flag)
  names(doubleRankedList)[which(!hits)] <- paste0(names(doubleRankedList)[which(!hits)], 1:sum(!hits))
  geneSet <- c(geneSetUp, geneSetDown)
  ans <- mwwGST(doubleRankedList, geneSet, moreDetails = moreDetails, alternative = "two.sided", verbose = verbose)
  ans$call <- ccall
  if(moreDetails) ans$doubleRankedList <- doubleRankedList
  ans$geneSetUp <- geneSetUp
  ans$geneSetDown <- geneSetDown
  class(ans) <- "mwwExtGST"
  invisible(ans)  
}
