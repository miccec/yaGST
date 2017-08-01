eeMWW <-
function(ddata, minoritySet, runs = 1000) {
  m <- length(minoritySet)
  doubled_m <- 2 * m
  majoritySet <- setdiff(rownames(ddata), minoritySet)
  
  a_foreach  <- foreach(k = 1:runs) %dopar% {
    sampled_majoritySet <- sample(majoritySet, doubled_m, replace = FALSE)
    tmp <- apply(ddata, 2, function(x) wilcox.test(x[minoritySet], x[sampled_majoritySet])$statistic)
    return(tmp)
  }
  ans <- a_foreach[[1]]
  for(k in 2:length(a_foreach)) ans <- ans + a_foreach[[k]]  
  ans <- sort(ans/runs/m/doubled_m, decreasing = TRUE)
  return(ans)
}
