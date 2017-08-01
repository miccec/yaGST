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
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
library(mwwGSEA)
data(rankedList)
head(rankedList, 10)
tail(rankedList, 10)
fivenum(rankedList)
}
\keyword{datasets}
