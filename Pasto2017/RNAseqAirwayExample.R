
##Usaremos 8 muestras del paquete airway, el cual hace parte del articulo
##[Himes et al](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4057123/): 
##"RNA-seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates 
## Cytokine Function in Airway Smooth Muscle Cells".

# First, make variables for the different BAM files and GTF file. Use the `sample.table`
# to contruct the BAM file vector, so that the count matrix will be in the same order as the `sample.table`.

library(airway)
data(airway)
airway
colData(airway)
rowRanges(airway)

##The counts matrix is stored in `assay` of a *SummarizedExperiment*.
head(assay(airway))
library(rafalib)
mypar()
## Note that, on the un-transformed scale, the high count genes have high variance.
## Cone shape; This is a general property of counts generated from sampling processes, 
## that the variance typically increases with the expected value.
plot(assay(airway)[,1:2], cex=.1,main="counts")

## Working with DESeq; The *DESeqDataSet* object 
## We need to specify an experimental *design* here, for later use in differential analysis. 
## The design starts with the tilde symbol `~`, which means, model the counts (log2 scale) 
## using the following formula. Following the tilde, the variables are columns of the `colData`, 
## and the `+` indicates that for differential expression analysis we want to compare levels of `dex` 
## while controlling for the `cell` differences.
library(DESeq2)
dds <- DESeqDataSet(airway, design= ~ cell + dex)

###############################################################################
### Optional if we wanto to run from bam files and get the count matrix ##
library(airway)
dir <- system.file("extdata", package="airway", mustWork=TRUE)
csv.file <- file.path(dir, "sample_table.csv")
sample.table <- read.csv(csv.file, row.names=1)
bam.files <- file.path(dir, paste0(sample.table$Run, "_subset.bam"))
gtf.file <- file.path(dir, "Homo_sapiens.GRCh37.75_subset.gtf")
## construct the count matrix
library(Rsubread)
fc <- featureCounts(bam.files, annot.ext=gtf.file,
                    isGTFAnnotationFile=TRUE, 
                    isPaired=TRUE)
names(fc)
unname(fc$counts) # hide the colnames
##  The *DESeqDataSet* object if we have a count matrix
colnames(fc$counts) <- rownames(sample.table)
dds.fc <- DESeqDataSetFromMatrix(fc$counts, 
                                 colData=sample.table, 
                                 design=~ cell + dex)

#############################################################################

###* Normalization for sequencing depth
## The following estimates size factors to account for differences in sequencing depth, 
## and is only necessary to make the `log.norm.counts` object below.

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
colSums(counts(dds))
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

#################################################################################
### Size factors manual:
## Size factors are calculated by the median ratio of samples to a pseudo-sample 
## (the geometric mean of all samples). In other words, for each sample, 
## we take the exponent of the median of the log ratios in this histogram.

loggeomeans <- rowMeans(log(counts(dds)))
hist(log(counts(dds)[,1]) - loggeomeans, 
     col="grey", main="", xlab="", breaks=40)

# The size factor for the first sample:
exp(median((log(counts(dds)[,1]) - loggeomeans)[is.finite(loggeomeans)]))
sizeFactors(dds)[1]
#########################################################################

## Make a matrix of log normalized counts (plus a pseudocount):
log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
## or
log.norm <- normTransform(dds)

## Examine the log counts and the log normalized counts (plus a pseudocount).

rs <- rowSums(counts(dds))
mypar(1,2)
boxplot(log2(counts(dds)[rs > 0,] + 1)) # not normalized
boxplot(log.norm.counts[rs > 0,]) # normalized

## Scatterplot: Note the fanning out of the points in the lower left corner, for points less than 2^5
mypar(1,3)
plot(log.norm.counts[,1:2], cex=.1,main="log-norm counts")

### Stabilizing count variance
## Now we will use a more sophisticated transformation, which is similar to the variance stablizing normalization method
## It uses the variance model for count data to shrink together the log-transformed counts for genes with very low counts.
## For genes with medium and high counts, the `rlog` is very close to `log2`
rld <- rlog(dds)
plot(assay(rld)[,1:2], cex=.1,main="rlog-norm counts")
## Another transformation for stabilizing variance in the *DESeq2* package is `varianceStabilizingTransformation`. 
## These two tranformations are similar, the *rlog* might perform a bit better when the size factors vary widely, 
## and the *varianceStabilizingTransformation* is much faster when there are many samples.
vsd <- varianceStabilizingTransformation(dds)
plot(assay(vsd)[,1:2], cex=.1,main="vsdlog-norm counts")

## We can examine the standard deviation of rows over the mean. 
## Note that the genes with high variance for the *log* come from the genes with lowest mean. 
## If these genes were included in a distance calculation, 
## the high variance at the low count range might overwhelm the signal at the higher count range.
library(vsn)
meanSdPlot(log.norm.counts, ranks=FALSE) 
##For the rlog:
meanSdPlot(assay(rld), ranks=FALSE)
##For the VST:
meanSdPlot(assay(vsd), ranks=FALSE)

### Visulaize difference and similarities among samples
## The principal components (PCA) plot is a useful diagnostic for examining relationships between samples:
plotPCA(log.norm, intgroup="dex")
#Using the rlog:
plotPCA(rld, intgroup="dex")
#Using the VST:
plotPCA(vsd, intgroup="dex")
## clustering
mypar(1,2)
plot(hclust(dist(t(log.norm.counts))), labels=colData(dds)$dex)
plot(hclust(dist(t(assay(rld)))), labels=colData(rld)$dex)

######### ****** DE ******* #############

## Differential gene expression

### Modeling raw counts with normalization

## We will use an overdispersed Poisson distribution -- called the negative binomial -- to model the *raw counts* 
## in the count matrix. The model will include the *size factors* into account to adjust for sequencing depth and the *dispersion parameter*.

mypar(3,1)
n <- 10000
brks <- 0:300
hist(rpois(n,100),main="",xlab="",breaks=brks,col="black")
hist(rpois(n,1000)/10,main="",xlab="",breaks=brks,col="black")
hist(rpois(n,10)*10,main="",xlab="",breaks=brks,col="black")

# So, when we scale a raw count, we break the implicit link between the mean and the variance. 
# This is not necessarily a problem, if we have 100s of samples over which to observe within-group variance, 
# however RNA-seq samples can often have only 3 samples per group, in which case, 
# we can get a benefit of information from using raw counts, and incorporating normalization factors 
# on the right side of the equation above.

### Counts across biological replicates and over-dispersion

## For the negative binomial, the variance parameter is called *disperison*, 
## and it links the mean value with the expected variance. 
## The reason we see more dispersion than in a Poisson is mostly due to changes in the 
## proportions of genes across biological replicates -- which we would expect 
## due to natural differences in gene expression. 

mypar(3,1)
n <- 10000
brks <- 0:400
hist(rpois(n,lambda=100),
     main="Poisson / NB, disp=0",xlab="",breaks=brks,col="black")
hist(rnbinom(n,mu=100,size=1/.01),
     main="NB, disp = 0.01",xlab="",breaks=brks,col="black")
hist(rnbinom(n,mu=100,size=1/.1),
     main="NB, disp = 0.1",xlab="",breaks=brks,col="black")

## The square root of the dispersion is the coefficient of variation -- SD/mean -- 
## after subtracting the variance we expect due to Poisson sampling.

disp <- 0.5
mu <- 100
v <- mu + disp * mu^2
sqrt(v)/mu
sqrt(v - mu)/mu
sqrt(disp)

### Diff

levels(dds$dex)
dds$dex <- relevel(dds$dex, "untrt")
levels(dds$dex)

dds <- DESeq(dds)
res <- results(dds)
### Examining results tables
head(res)
table(res$padj < 0.1)
#A summary of the results can be generated:
summary(res)
#For testing at a different threshold, we provide the `alpha` to *results*, 
# so that the mean filtering is optimal for our new FDR threshold.
res2 <- results(dds, alpha=0.05)
table(res2$padj < 0.05)

##############################
###* Visualizing results *###
############################
# Los MA-plot proporcionan una vista global de los genes diferenciales, 
# en el eje "y" el log2 fold change y sobr el eje "x" la media de los conteos normalizados
mypar()
plotMA(res, ylim=c(-4,4))
plotMA(res2, ylim=c(-4,4))

# Podemos visualizar aquellos que tienen un fold change mas o menos de el doble de la expresión
# que los controles
res.thr <- results(dds, lfcThreshold=1)
plotMA(res.thr, ylim=c(-4,4))

# Histograma del p-value
hist(res$pvalue[res$baseMean > 1], 
     col="grey", border="white", xlab="", ylab="", main="")

# Una tabla de resultados ordenada de acuerdo al pvalue
resSort <- res[order(res$padj),]
head(resSort)

# Examinar los conteos para el genes con mejor adj.p-value
plotCounts(dds, gene=which.min(res$padj), intgroup="dex")

# plot de conteos para el gen con mejor adj.p-value  con ggplot:
library(ggplot2)
data <- plotCounts(dds, gene=which.min(res$padj), intgroup=c("dex","cell"), returnData=TRUE)
ggplot(data, aes(x=dex, y=count, col=cell)) +
  geom_point(position=position_jitter(width=.1,height=0)) +
  scale_y_log10()

ggplot(data, aes(x=dex, y=count, col=cell, group=cell)) +
  geom_point() + geom_line() + scale_y_log10() 

library(pheatmap)
topgenes <- head(rownames(resSort),20)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("dex","cell")])
pheatmap(mat, annotation_col=df)


library(sva)
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ dex, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=2)

#Do the surrogate variables capture the cell difference?

plot(svseq$sv[,1], svseq$sv[,2], col=dds$cell, pch=16)

dds.sva <- dds
dds.sva$SV1 <- svseq$sv[,1]
dds.sva$SV2 <- svseq$sv[,2]
design(dds.sva) <- ~ SV1 + SV2 + dex
dds.sva <- DESeq(dds.sva)
res.sva <- results(dds.sva)

hist(res$pvalue[res$baseMean > 1], 
     col="grey", border="white", xlab="", ylab="", main="")

hist(res.sva$pvalue[res.sva$baseMean > 1], 
     col="grey", border="white", xlab="", ylab="", main="")

