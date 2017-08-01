print.mwwExtGST <-
function (obj) {
  ans <- list()
  if (length(obj) < 2) {
    message("Nothing to print.")
    return(invisible(NULL))
  }
  
  print(unclass(obj$call))
  message(paste0("length of the original positive gene-set: ", length(obj$geneSetUp)))
  message(paste0("length of the original negative gene-set: ", length(obj$geneSetDown)))
  message(paste0("length of the pooled gene-set: ", obj$originalGeneSetCount))
  message(paste0("length of the pooled gene-set matched with the genes in the gene-list: ", obj$actualGeneSetCount))
  message(paste0("length of the doubled gene-list: ", obj$lengthOfRankedList))
  message(paste0("NES = ", round(obj$nes, 4), ", prob.unbalance = ", round(obj$pu, 4), ", log2(prob.unbalance) = ", round(obj$log.pu, 4)))
  message(paste0("p-value = ", format(obj$p.value, digits = 4), ", for the alternative hypothesis: ", obj$alternative))
  invisible(NULL)
}
