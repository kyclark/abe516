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
library("biomaRt")
library("pheatmap")
library("RColorBrewer")
library("reshape2")
library("dplyr")
library("genefilter")

wd = "~/work/abe516-project2"
fig.dir = file.path(wd, "figures")
if (!dir.exists(fig.dir)) {
  dir.create(fig.dir)
}

filenames = dir(path = file.path(wd, "stage4_bamfiles"), pattern = "bam", full.names = T)
printf("Found %s BAM files\n", length(filenames))
stage4.samples = read.csv(file.path(wd, "stage4_samples.csv"), header = F)
colnames(stage4.samples) = c("SampleName", "trt")

gtf.file = file.path(wd, "sorghum.gtf")
txdb = makeTxDbFromGFF(gtf.file, format = "gtf", circ_seqs = character())
ebg = exonsBy(txdb, by = "gene")

# this function helps to prepare data if we have more filenames and sample tables 
dataprepare = function(filenames, sample.table) {
  bam.files = BamFileList(filenames)
  se = summarizeOverlaps(features = ebg, reads = bam.files, mode = "Union", singleEnd = TRUE)
  colData(se) <- DataFrame(sample.table)
  return(se)
}

se = dataprepare(filenames = filenames, sample.table = stage4.samples)
se
dim(se)
assayNames(se)
head(assay(se))
colSums(assay(se))

se$trt
countdata = assay(se)
head(countdata)
coldata = colData(se)

dds = DESeqDataSet(se, design = ~trt)
nrow(dds)
dds = dds[rowSums(counts(dds)) > 1, ]  ### pre-filtering
nrow(dds)

# PCA 
rld = rlog(dds, blind = FALSE)  ## transform data to log2 scale
pcdata = plotPCA(rld, intgroup = "trt", returnData = TRUE)
percentVar = round(100 * attr(pcdata, "percentVar"))

pca.plot = ggplot(pcdata, aes(PC1, PC2, color = trt, shape = trt)) + 
  ggtitle("Stage4") +
  geom_point(size = 3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

ggsave(file = file.path(fig.dir, "pca.png"), width = 5, height = 5, plot = pca.plot)

# Heatmap
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$trt
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         main = "Sample Distances",
         filename = file.path(fig.dir, "heatmap.png"),
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# MDS plot
mds = as.data.frame(colData(rld)) %>% cbind(cmdscale(sampleDistMatrix))
mds.plot = ggplot(mds, aes(x = `1`, y = `2`, color = trt, shape = trt)) + 
  ggtitle("MDS Plot") + geom_point(size = 3) + coord_fixed() + theme_light()
ggsave(file = file.path(fig.dir, "mds.png"), plot = mds.plot)

# Dendrogram
dist.matrix = as.dist(1 - sampleDistMatrix)
fit = hclust(dist.matrix, method = "ward.D2") 
dg = ggdendro::ggdendrogram(fit, rotate=T) + ggtitle("Dendrogram")
ggsave(file = file.path(fig.dir, "dendrogram.png"), width = 5, height = 5, plot = dg)

# Gene clustering
topVarGenes = head(order(rowVars(assay(rld)), decreasing = TRUE), 25)
mat  = assay(rld)[ topVarGenes, ]
mat  = mat - rowMeans(mat)
colnames(mat) = rld$trt
pheatmap(mat, main="Gene Clusters", filename = file.path(fig.dir, "gene_clusters.png"))

# Plot GO functions
go = read.csv(file = file.path(wd, "go.csv"), header = T)
go.terms = unlist(strsplit(as.character(go$GO.term.name), ', ', fixed = T))
go.summary = as.data.frame(table(go.terms))
colnames(go.summary) = c("func", "freq")
go.summary$func = factor(go.summary$func, 
                         levels = go.summary$func[order(go.summary$freq)])
go.plot = ggplot(data = go.summary, aes(x=func, y=freq)) + xlab(element_blank()) +
  geom_bar(stat="identity", fill="#56B4E9", colour="black") + theme_light() + coord_flip()
ggsave(file = file.path(fig.dir, "go.png"), width = 5, height = 6, plot = go.plot)

dds = DESeq(dds)
res = results(dds)
mcols(res, use.names = TRUE)
summary(res)

resLFC1 = results(dds, lfcThreshold = 1)
table(resLFC1$padj < 0.1)

sum(res$padj < 0.05, na.rm = TRUE)

resSig = subset(res, padj < 0.05)
head(resSig[order(resSig$log2FoldChange), ])

resSig = resSig[order(resSig$log2FoldChange, decreasing = TRUE), ]
write.csv(resSig, file = file.path(wd, "gene_list.csv"))

# Create plots for the top 10 genes
for (i in 1:10) {
  gene.name = rownames(resSig)[i]
  printf("%s: Creating plot for '%s'\n", i, gene.name)
  d = plotCounts(dds, 
                 gene = gene.name, 
                 intgroup = c("trt"), 
                 returnData=TRUE)
  p = ggplot(d, aes(x = trt, y = count, color = trt, shape = trt)) + 
    ggtitle(gene.name) + theme_light() +
    geom_point(position = position_jitter(w = 0.1, h = 0))
  filename = file.path(fig.dir, 
                       paste0(paste("top", sprintf("%02d", i), gene.name, sep="_"), ".png"))
  ggsave(file = filename, width = 5, height = 5, plot = p)
}

# MAplot
png(filename = file.path(fig.dir, "ma_plot.png"))
plotMA(res, ylim = c(-3, 3), main = "MA Plot")
invisible(dev.off())

# plotmalabel
plotMA(resLFC1, ylim = c(-3, 3))
topGene = rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
  printf("%s\n", baseMean)
  printf("%s\n", log2FoldChange)
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, topGene, pos = 2, col = "dodgerblue")
})