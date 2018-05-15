library("R.utils")
library("AnnotationHub")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicRanges")
library("Rsubread")
library("limma")
library("edgeR")
library("DESeq2")
library("GenomicAlignments")
library("BiocParallel")
library("AnnotationDbi")
library("ggplot2")
library("Gviz")
library("org.Hs.eg.db")
library("biomaRt")
library("ReportingTools")

wd = "~/work/abe516-project2"
fig.dir = file.path(wd, "figures")
if (!dir.exists(fig.dir)) {
  dir.create(fig.dir)
}
#filenames = dir(path = file.path(wd, "bamfiles"), full.names = T)
filenames = dir(path = file.path(wd, "bamfiles"), pattern = "bam", full.names = T)
printf("Found %s BAM files\n", length(filenames))
sample.table = read.csv(file.path(wd, "sample_table.csv"), header = F)
colnames(sample.table) = c("SampleName", "trt")

fn_S4_SS = filenames[grep("stage_4_SS",sample.table$trt)]#subset stage 4 SS
fn_S4_PS = filenames[grep("stage_4_PS",sample.table$trt)]#subset stage 4 PS
st_S4_SS = sample.table[grep("stage_4_SS",sample.table$trt),]# subset stage 4 SS
st_S4_PS = sample.table[grep("stage_4_PS",sample.table$trt),]# subset stage 4 PS



gtf.file = file.path(wd, "sorghum.gtf")
txdb = makeTxDbFromGFF(gtf.file, format = "gtf", circ_seqs = character())
ebg = exonsBy(txdb, by = "gene")

# this function helps to prepare data if we have more filenames and sample tables 
dataprepare = function(filenames,sample.table){
  bam.files = BamFileList(filenames)
  se = summarizeOverlaps(features = ebg, 
                         reads = bam.files, 
                         mode = "Union", 
                         singleEnd = TRUE)
  colData(se) <- DataFrame(sample.table)
  return(se)
}

se = dataprepare(filenames = fn_S4_SS, sample.table = st_S4_SS)
se
dim(se)
assayNames(se)
head(assay(se))
colSums(assay(se))
  
# -------------------------------------------------------

#se$trt = relevel(se$trt, "control")
se$trt
countdata = assay(se)
head(countdata)
coldata = colData(se)

# -------------------------------------------------------

dds = DESeqDataSet(se, design = ~trt)
nrow(dds)
dds = dds[rowSums(counts(dds)) > 1, ]  ### pre-filtering
nrow(dds)

# -------------------------------------------------------
rld = rlog(dds, blind = FALSE)  ## transform data to log2 scale
pcdata = plotPCA(rld, intgroup = "trt", returnData = TRUE)
percentVar = round(100 * attr(pcdata, "percentVar"))

pca.plot = ggplot(pcdata, aes(PC1, PC2, color = trt, shape = trt)) + 
  geom_point(size = 3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

ggsave(file = file.path(fig.dir, "pca.png"), width = 5, height = 5, plot = pca.plot)

# -------------------------------------------------------
dds = DESeq(dds)  ### input data: raw count data
res = results(dds)
mcols(res, use.names = TRUE)

# -------------------------------------------------------
summary(res)

resLFC1 = results(dds, lfcThreshold = 1)
table(resLFC1$padj < 0.1)

sum(res$padj < 0.05, na.rm = TRUE)

resSig = subset(res, padj < 0.05)
head(resSig[order(resSig$log2FoldChange), ])

resSig = resSig[order(resSig$log2FoldChange, decreasing = TRUE), ]
write.csv(resSig, file = file.path(wd, "gene_list.csv"))

# plotcounts
topGene = rownames(res)[which.min(res$padj)]
png(filename = file.path(fig.dir, "top_gene.png"))
plotCounts(dds, gene = topGene, intgroup = c("trt"))
invisible(dev.off())

# MAplot
png(filename = file.path(fig.dir, "ma_plot.png"))
plotMA(res, ylim = c(-3, 3))
invisible(dev.off())


# plotmalfc plotMA(resLFC1, ylim=c(-3,3))

# plotmalabel
plotMA(resLFC1, ylim = c(-3, 3))
topGene = rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
    points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
    text(baseMean, log2FoldChange, topGene, pos = 2, col = "dodgerblue")
})

# What do we use for sorghum?
#columns(org.Hs.eg.db)
#columns(biomaRt)

ah <- AnnotationHub()
query(ah, "Sorghum") #check if sorghum in the dataset
# 2 resutls: AH10585 AH59058
ath <- ah[['AH59058']]txdb #retrieve dataset AH59058

keytypes(ath) #tell us which can be set as keytype

res=res[-seq(1:10),]
res$symbol = mapIds(ath, keys = row.names(res), column = "SYMBOL", keytype = "SORGHUM_BICOLOR", multiVals = "first")
res$entrez = mapIds(ath, keys = row.names(res), column = "ENTREZID", keytype = "GENENAME", multiVals = "first")

resOrdered = res[order(res$padj), ]
head(resOrdered)

(resGR = results(dds, lfcThreshold = 1, format = "GRanges"))

resGR$symbol = mapIds(ath, names(resGR)[1:5], "SYMBOL", "ENTREZID")

window = resGR[topGene] + 1e+06
strand(window) = "*"
resGRsub = resGR[resGR %over% window]
naOrDup = is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group = ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

sig = factor(ifelse(resGRsub$padj < 0.1 & !is.na(resGRsub$padj), "sig", "notsig"))

options(ucscChromosomeNames = FALSE)
g = GenomeAxisTrack()
a = AnnotationTrack(resGRsub, name = "gene ranges", feature = sig)
d = DataTrack(resGRsub, data = "log2FoldChange", baseline = 0, type = "h", name = "log2 fold change", 
    strand = "+")
plotTracks(list(g, d, a), groupAnnotation = "group", notsig = "grey", sig = "hotpink")