\name{rankedList}
\alias{rankedList}
\docType{data}
\title{An example of pre-ranked list.}
%%\description{An example of pre-ranked list.}
\usage{data("rankedList")}
\format{
  The format is:
 Named num [1:17814] 4.82 4.33 4.25 4.18 4.09 ...
 - attr(*, "names")= chr [1:17814] "1.48043767313367" "1.37586701949352" "0.212142344466495" "-0.0291897442883637" ...
}
\details{This ordered list of genes from the comparison between G-CIMP-Low versus G-CIMP-High in the GBM. 
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  Ceccarelli at al - \emph{Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma.} - CELL, Volume 164, Issue 3, p550–563, 28 January 2016)
}


\examples{
library(yaGST)
data(rankedList)
head(rankedList, 10)
tail(rankedList, 10)
fivenum(rankedList)
}
\keyword{datasets}
