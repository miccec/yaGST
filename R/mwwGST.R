mwwGST <-
function (rankedList, geneSet, minLenGeneSet = 5, alternative = "greater", 
                     moreDetails = FALSE, verbose = TRUE) 
{
  ccall <- match.call()
  gs <- which(geneSet %in% names(rankedList))
  if (length(gs) < minLenGeneSet) {
    if (verbose) 
      message("The number of elements of the gene set is less than the minimum allowed.")
    ans <- list()
    class(ans) <- "mwwGST"
    return(ans)
  }
  gs <- geneSet[gs]
  outside_gs <- setdiff(names(rankedList), gs)
  ans <- list()
  ans$call <- ccall
  ans$alternative <- alternative
  ans$originalGeneSetCount <- length(geneSet)
  ans$actualGeneSetCount <- length(gs)
  ans$actualGeneSet <- gs
  ans$rankedList <- rankedList
  ans$lengthOfRankedList <- length(rankedList)
  tmp <- wilcox.test(rankedList[gs], rankedList[outside_gs], alternative = alternative)
  ans$statistic <- tmp$statistic; names(ans$statistic) <- NULL
  ans$nes <- tmp$statistic/length(gs)/length(outside_gs); names(ans$nes) <- NULL
  ans$pu <- ans$nes/(1 - ans$nes)
  ans$log.pu <- log2(ans$pu)
  ans$p.value <- tmp$p.value
  if (!moreDetails) 
    ans$rankedList <- NULL
  class(ans) <- "mwwGST"
  invisible(ans)
}
