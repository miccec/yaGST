gmt2GO <-
function(what) {
  tmp <- function(x) {    
    GO <- readLines(x)
    GO <- strsplit(GO, "\t")
    nnames <- lapply(GO, function(x) x[1])
    link <- lapply(GO, function(x) x[2])
    GO <- lapply(GO, function(x) x[3:length(x)])
    for(k in 1:length(GO)) {
      attr(GO[[k]], "link") <- link[k]
      attr(GO[[k]], "ontology") <- x
    }
    names(GO) <- nnames
    return(GO)
  }
  GO <- tmp(what[1])
  
  if(length(what) > 1)  
    for(k in 2:length(what)) GO <- c(GO, tmp(what[k]))
  
  invisible(GO)
}
