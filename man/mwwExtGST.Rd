\name{mwwExtGST}
\alias{mwwExtGST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Competitive Regulon Test}
\description{Run a competitive test to highlight whether a regulon, with postive and negative gene associated with a transcription factor, is highly ranked in a sequence of gene values.}
\usage{
mwwExtGST(rankedList, geneSetUp, geneSetDown, minLenGeneSet = 15, moreDetails = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{rankedList}{numeric vector of data values where the names are the genes names}
  \item{geneSetUp}{a character list of genes having a positive association with the transcription factor.}
  \item{geneSetDown}{a character list of genes having a negative association with the transcription factor.}
  \item{minLenGeneSet}{minimum dimension of the pooled geneSet}
  \item{moreDetails}{a logical indicating whether the output includes the rankedList (necessary to plot the enrichment)}
  \item{verbose}{a logical indicating to suppress or not the messages; it's TRUE by default.}
}
\details{The rankedList has to be a named sequence of values where the genes associated with the phenotype are positvive values, while those not associated are negative. This is necessary because the \verb{doubledRankedList} is set as \verb{c(rankedList, -rankedList)}.}

\value{
\item{call}{a character string of the call of the function.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{originalGeneSetCount}{the length of the pooled positive and negative gene-sets.}
\item{geneSetUp}{the same character list given in input.}
\item{geneSetDown}{the same character list given in input.}
\item{actualGeneSet}{the list of pooled positive and negative genes as comes from the intesection between the gene-set and the ranked-list.}
\item{actualGeneSetCount}{the length of the \verb{actualGeneSet}.}
\item{doubleRankedList}{the doubled ranked-list given in input; this slot is filled whether moreDetails is TRUE. Seeq details}
\item{lengthOfRankedList}{the length of the ranked-list given in input}
\item{statistic}{the value of the Mann-Whitney-Wilcox test statistic.}
\item{nes}{the value of the normalized enrichment score.}
\item{pu}{is the probability unbalance, i.e. the ratio of nes to 1-nes.}
\item{log.pu}{the log2 transformation of the pu.}
\item{p.value}{the p-value for the test.}
}

\references{
  Lim, W. K., Lyashenko, E. and Califano, A. - \emph{Master regulators used as breast cancer metastasis classifier.} - Pac Symp Biocomput, 504-515 (2009))
}
\author{Stefano M. Pagnotta}
\note{This function adatpts the enrichment analysis methodology from Lim at el. (2009) to the mwwGST function. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{rankedList}}, \code{\link{mwwGST}}}

\examples{
  library(yaGST)
  data("rankedList")
  positive_gs <- sample(head(names(rankedList), 10000), 200)
  negative_gs <- sample(tail(names(rankedList), 10000), 200)
  ans <- mwwExtGST(rankedList, positive_gs, negative_gs, moreDetails = TRUE)
  ans
  plot(ans)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
