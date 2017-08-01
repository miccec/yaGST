\name{eeMWW}
\alias{eeMWW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{eeMWW}
\description{This function implements the Easy Ensemble, together with the Mann-Witney-Wilcox test, to detect the genes associated with few samples (minority set) being a subset of a larger collection of samples (majority set).
}
\usage{
eeMWW(ddata, minoritySet, runs = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ddata}{a matrix where the samples are by rows and the features are in the columns.}
  \item{minoritySet}{a character vector of the minority set matching some row names of ddata.}
  \item{runs}{number of resampling.}
}
\details{The EasyEnsemble (EE) resampling scheme is an
Undersampling technique aimed to compare few samples (minority set), carrying some phenotype, to a larger collection of samples (majority set) unrelated with the phenotype. We implement the EE with the Mann-Whitney-Wilcoxon test (MWW) to compare the minority set, of dimension m, with a randomly selected collection of 2*m samples from the majority set. 
}
\value{a named vector of real values.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Xu-Ying Liu, Jianxin Wu, and Zhi-Hua Zhou - \emph{Exploratory Undersampling for Class-Imbalance Learning} - IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS?PART B: CYBERNETICS, VOL. 39, NO. 2, APRIL 2009}

\author{Stefano M. Pagnotta}
\note{We suggest running the function in a parallel setup.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{mwwGST}}}

\examples{require(yaGST)
nr <- 100; nc <- 1000
# generate a data-matrix with nr samples, and nc features
exprData <- matrix(rpois(nc * nr, 10), nrow = nr, ncol = nc)
colnames(exprData) <- paste0("feat", 1:nc)
rownames(exprData) <- paste0("sam", 1:nr)

# increase the first 3 samples (minority set) of 10\% of the original intensity 
# of the first 30 features (later the gene-set)
exprData[1, 1:30] <- exprData[1, 1:30]* runif(30, min = 1, max = 1.10)
exprData[2, 1:30] <- exprData[1, 1:30]* runif(30, min = 1, max = 1.10)
exprData[3, 1:30] <- exprData[1, 1:30]* runif(30, min = 1, max = 1.10)
samples_of_interest <- rownames(exprData)[1:3] # minority set

# running in parallel
library(doParallel)
# adjust the number of CPUs as needed
cl <- makePSOCKcluster(3)
clusterApply(cl, floor(runif(length(cl), max = 10000000)), set.seed)
registerDoParallel(cl)
ans_eeMWW <- eeMWW(exprData, samples_of_interest)
stopCluster(cl)

# set the gene-set and run the enrichment analysis
geneSet <- colnames(exprData)[1:30]
(tmp <- mwwGST(ans_eeMWW, geneSet))
plot(tmp, rankedList = ans_eeMWW)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
