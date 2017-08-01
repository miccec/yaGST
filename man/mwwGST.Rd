\name{mwwGST}
\alias{mwwGST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Competitive Gene Set Test}
\description{Run a competitive test to highlight whether a gene set is highly ranked in a sequence of gene values on the genes outside the gene-set.}
\usage{
mwwGST(rankedList, geneSet, minLenGeneSet = 5, alternative = "greater", moreDetails = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rankedList}{numeric vector of data values where the names are the genes names}
  \item{geneSet}{a character list of genes}
  \item{minLenGeneSet}{minimum dimension of the geneSet}
  \item{alternative}{a character string specifying the alternative hypothesis ("two.sided", "less", "greater").}
  \item{moreDetails}{a logical indicating whether the output includes the rankedList (necessary to plot the enrichment)}
  \item{verbose}{a logical indicating to suppress or not the messages; it's TRUE by default.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\item{call}{a character string of the call of the function.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{originalGeneSetCount}{the length of the gene-set given in input}
\item{actualGeneSet}{the list of genes as comes from the intesection between the gene-set and the ranked-list.}
\item{actualGeneSetCount}{the length of the actualGeneSet}
\item{rankedList}{the ranked-list given in input; this slot is filled whether moreDetails is TRUE}
\item{lengthOfRankedList}{the length of the ranked-list given in input}
\item{statistic}{the value of the Mann-Whitney-Wilcox test statistic.}
\item{nes}{the value of the normalized enrichment score.}
\item{pu}{is the probability unbalance, i.e. the ratio of nes to 1-nes.}
\item{log.pu}{the log2 transformation of the pu.}
\item{p.value}{the p-value for the test.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{Stefano M. Pagnotta}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{rankedList}}, \code{\link{mwwExtGSEA}}}
\examples{
  library(yaGST)
  data("rankedList")
  # generate a random data set of dimension 100
  geneSet <- sample(head(names(rankedList), 5000), 100)
  ans <- mwwGST(rankedList, geneSet, moreDetails = TRUE)
  ans
  plot(ans)
  
  # generate a second gene set
  geneSet <- sample(tail(names(rankedList), 5000), 100)
  ans <- mwwGST(rankedList, geneSet, moreDetails = TRUE)
  plot(ans)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
