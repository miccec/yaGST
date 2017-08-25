# yaGST

This is a collection of wrappers to the Wilcoxon test to run competitive gene set tests and comparison between imbalanced datasets.


## Installation

```{r}
library(devtools)
install_github("miccec/yaGST")
```
## Running a gene-set enrichment analysis

```{r}
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

```
## How to run an enrichment analysis on a collection of gene-sets.
You may need to download the collection of gene-sets, for example from the Molecular Signatures Database at [the Broad Institute](http://software.broadinstitute.org/gsea/downloads.jsp)

- [h.all.v6.0.symbols.gmt](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.0/h.all.v6.0.symbols.gmt)

In case the output of the analysis has to be saved in an .xls file, you need to install the [WriteXLS](https://cran.r-project.org/web/packages/WriteXLS/index.html) package.

load the GO collections

```{r}
library(yaGST)
GO <- gmt2GO("~/pathToYourFiles/h.all.v6.0.symbols.gmt")
length(GO)
```
load the ranked list
```{r}
data("rankedList")
head(rankedList)
tail(rankedList)
```
run the mwwGST
```{r}
ans <- lapply(GO, function(x) mwwGST(rankedList, x, minLenGeneSet = 15, verbose = FALSE))
```
Now, arrange the results in a table

```{r}
actualGeneSetSize <- unlist(sapply(ans, function(x) x$actualGeneSetCount))
originalGeneSetSize <- unlist(sapply(ans, function(x) x$originalGeneSetCount))
pValue <- unlist(sapply(ans, function(x) x$p.value))
qValue <- p.adjust(pValue)
nes <- unlist(sapply(ans, function(x) x$nes))
pu <- unlist(sapply(ans, function(x) x$pu))
log2.pu <- unlist(sapply(ans, function(x) x$log.pu))
result_table <- cbind(originalGeneSetSize, actualGeneSetSize,  
  round(nes, 4), round(pu, 4), round(log2.pu, 4), round(pValue, 6), round(qValue, 6))
colnames(result_table) <- c("originalGeneSetSize", "actualGeneSetSize", "NES", "pu", "log2_pu", "pValue", "qValue")
result_table <- result_table[order(result_table[, "NES"], decreasing = TRUE),]
```


If needed, the table can be stored in an .xls file as follows 

```{r, eval = FALSE}
library(WriteXLS)
gst_results <- as.data.frame(result_table)
WriteXLS("gst_results", ExcelFileName = "~/pathToYourFiles/aNameForTheAnalysis.xls", row.names = TRUE)
```

filter and plot the results
```{r}
qValue_treshold <- 0.05
NES_treshold <- 0.6
selected <- which((result_table[, "NES"] > NES_treshold) & (result_table[, "qValue"] < qValue_treshold))
result_table[selected,]
first_geneSet <- names(selected)[1]
plot(ans[[first_geneSet]], rankedList = rankedList, main = first_geneSet)
```

## Running an extended enrichment analysis on a regulon (with positive and negative targets)
```{r}
data("rankedList")
positive_gs <- sample(head(names(rankedList), 10000), 200)
negative_gs <- sample(tail(names(rankedList), 10000), 200)
ans <- mwwExtGST(rankedList, positive_gs, negative_gs, moreDetails = TRUE)
ans
plot(ans)
```

## Running ee-MWW

```{r}
nr <- 100; nc <- 1000
# generate a data-matrix with nr samples, and nc features
exprData <- matrix(rpois(nc * nr, 100), nrow = nr, ncol = nc)
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
```
