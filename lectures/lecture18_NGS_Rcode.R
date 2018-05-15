# To run this case study, you need have R version of 3.2.1 or later

##
##########################################################################

# Aligning reads to genome using Bowtie2

#The first step for reads alignment is to download reference index.  Most reference index files can be downloaded from URL https://support.illumina.com/sequencing/sequencing_software/igenome.html
######################################################################

# load libraries
source("http://bioconductor.org/biocLite.R")

# Each of these commands tells Bioconductor to download and install a package
biocLite("GenomicRanges")
biocLite("GenomicFeatures")
biocLite("Gviz")
biocLite("Rsamtools")
biocLite("DESeq2")
biocLite("edgeR")  ### you may already have it. 
## biocLite("org.Mm.eg.db")     ### mouse sequence
biocLite("org.Hs.eg.db")   ### human sequence
biocLite("limma")  ### you may already have it. 
biocLite("Rsubread")
## biocLite("readGAlignmentsFromBam")
biocLite("GenomicAlignments")

library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(Rsubread)
library(limma)
library(edgeR)
library(DESeq2)
library(GenomicAlignments)
library(BiocParallel)
library("AnnotationDbi")
library("org.Hs.eg.db")
library("ggplot2")
library("Gviz")

wd = "~/work/abe516/rna-seq"
filenames <- dir(path = wd, pattern = "bam", full.names = T)
file.exists(filenames)
sampleTable <- read.csv(file.path(wd, "SraRunInfo.csv"), header = T)
## ------------------------------------------------------------------------

bamfiles <- BamFileList(filenames)

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
## transfer the file from HPC to your local computer :
## login to your HPC account: ssh yourNetId@hpc.arizona.edu
## cp /genome/iGenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf ~/temp/
## and then in your local computer (in the same directory as .bam files): 
## scp yourNetId@filexfer.hpc.arizona.edu:~/temp/genes.gtf .


gtffile <-dir(pattern="genes.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())
ebg <- exonsBy(txdb, by="gene")
## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
se <- summarizeOverlaps(features=ebg, 
                        reads=bamfiles,
                        mode="Union",
                        singleEnd=TRUE)
## ------------------------------------------------------------------------
se
dim(se)
assayNames(se)
head(assay(se))
colSums(assay(se))

## ------------------------------------------------------------------------
colData(se) <- DataFrame(sampleTable)
## ------------------------------------------------------------------------
se$trt <- relevel(se$trt, "control")
se$trt
countdata <- assay(se)
head(countdata)
coldata <- colData(se)
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

dds <- DESeqDataSet(se, design = ~  trt)
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1, ]  ### pre-filtering
nrow(dds)

## ----rlog----------------------------------------------------------------
rld <- rlog(dds, blind=FALSE)  ## transform data to log2 scale


pcdata <- plotPCA(rld, intgroup = "trt", returnData=TRUE)
percentVar <- round(100 * attr(pcdata, "percentVar"))

ggplot(pcdata, aes(PC1, PC2, color=trt, shape=trt)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

## ---------------------------------------------------------------
dds <- DESeq(dds)    ### input data: raw count data
res <- results(dds)
mcols(res, use.names=TRUE)

## ------------------------------------------------------------------------
summary(res)

## ------------------------------------------------------------------------
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

## ------------------------------------------------------------------------
sum(res$padj < 0.05, na.rm=TRUE)

## ------------------------------------------------------------------------
resSig <- subset(res, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])

## ------------------------------------------------------------------------
resSig=resSig[order(resSig$log2FoldChange, decreasing=TRUE), ]
write.csv(resSig, file="Myresult.csv")

## ----plotcounts----------------------------------------------------------
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("trt"))

## ----MAplot--------------------------------------------------------------
plotMA(res, ylim=c(-3,3))

## ----plotmalfc-----------------------------------------------------------
## plotMA(resLFC1, ylim=c(-3,3))

## ----plotmalabel---------------------------------------------------------
plotMA(resLFC1, ylim=c(-3,3))
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
columns(org.Hs.eg.db)

## ------------------------------------------------------------------------
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

## ------------------------------------------------------------------------
resOrdered <- res[order(res$padj),]
head(resOrdered)
## ------------------------------------------------------------------------
(resGR <- results(dds, lfcThreshold=1, format="GRanges"))

## ------------------------------------------------------------------------
resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")

## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

## ------------------------------------------------------------------------
sig <- factor(ifelse(resGRsub$padj < .1 & !is.na(resGRsub$padj),"sig","notsig"))

## ----gvizplot------------------------------------------------------------
options(ucscChromosomeNames=FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name="gene ranges", feature=sig)
d <- DataTrack(resGRsub, data="log2FoldChange", baseline=0,
               type="h", name="log2 fold change", strand="+")
plotTracks(list(g,d,a), groupAnnotation="group", notsig="grey", sig="hotpink")

################# The END ########################
