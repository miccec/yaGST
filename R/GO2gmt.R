GO2gmt <-
function(GO_, fileName) {
  geneSetName <- names(GO_)
  link  <- unlist(lapply(GO_, function(x) attr(x, "link")))
  geneSet <- lapply(GO_, function(x) paste0(x, collapse = '\t'))
  tmp <- paste(geneSetName, link, geneSet, sep = "\t")
  cat(tmp[1], file = fileName, append = FALSE, sep = "\n")
  for(k in 2:length(tmp)) cat(tmp[k], file = fileName, append = TRUE, sep = "\n")
  return(NULL)
}
