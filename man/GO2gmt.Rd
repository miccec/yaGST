\name{GO2gmt}
\alias{GO2gmt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a .gmt file from a list of gene set sequences.}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GO2gmt(GO_, fileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GO_}{a named vector list.
%%     ~~Describe \code{GO_} here~~
}
  \item{fileName}{a character string naming a file
%%     ~~Describe \code{fileName} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Stefano M. Pagnotta}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{gmt2GO}}, \code{\link{rankedList}}}

\examples{
library(mwwGSEA)
data("rankedList")

# create a collection of gene sets
GO <- vector("list", 2)
GO[[1]] <- sample(head(names(rankedList), 5000), 50)
# your reference lik for the gene set
attr(GO[[1]], "link") <- "http://www.enjoy_the_silence.dm" 

GO[[2]] <- sample(head(names(rankedList), 5000), 50)
attr(GO[[2]], "link") <- "http://www.imagine.jl"
names(GO) <- c("geneSet_1", "geneSet_2")
GO

# save the collection 
GO2gmt(GO, "~/my_GO_collection.gtm")
#########
# load a .gmt file
my_GO_collection <- gmt2GO("~/my_GO_collection.gtm")
summary(my_GO_collection)
head(my_GO_collection$geneSet_1)
attr(my_GO_collection[[1]], "link")
attr(my_GO_collection[[1]], "ontology")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
