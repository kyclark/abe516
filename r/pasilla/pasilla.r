#source("https://bioconductor.org/biocLite.R")
#biocLite("pasilla")
#biocLite("mixOmics")
#biocLite("DESeq2")

library("pasilla")
library("ggplot2")
library("reshape")
library("mixOmics")
library("RColorBrewer")
library("DESeq2")

datafile = system.file("extdata/pasilla_gene_counts.tsv", package = "pasilla")
rawCountTable = read.table(datafile, header = TRUE, row.names = 1)
condition = c("control", "control", "control", "control", "treated", "treated", "treated")
libType = c("single-end", "single-end", "paired-end", "paired-end", "single-end",
            "paired-end", "paired-end")
colnames(rawCountTable)[1:4] = paste0("control", 1:4)

# raw plot
ggplot(rawCountTable, aes(x = control1)) + geom_histogram(fill = "#525252", binwidth = 2000)

# log-normalized
pseudoCount = log2(rawCountTable + 1)
ggplot(pseudoCount, aes(x = control1)) + ylab(expression(log[2](count + 1))) +
  geom_histogram(colour = "white", fill = "#525252", binwidth = 0.6)

# boxplots
df = melt(pseudoCount, variable_name = "Samples")
df = data.frame(df, Condition = substr(df$Samples, 1, 7))
ggplot(df, aes(x = Samples, y = value, fill = Condition)) + geom_boxplot() + xlab("") +
  ylab(expression(log[2](count + 1))) + scale_fill_manual(values = c("#619CFF", "#F564E3"))

# histogram
ggplot(df, aes(x = value, colour = Samples, fill = Samples)) + ylim(c(0, 0.25)) +
  geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ Condition) +
  theme(legend.position = "top") + xlab(expression(log[2](count + 1)))

# MA plot
x = pseudoCount[, 1]
y = pseudoCount[, 2]

# M-values
M = x-y 

# A-values
A = (x + y)/2
df = data.frame(A, M)
ggplot(df, aes(x = A, y = M)) + geom_point(size = 1.5, alpha = 1/5) +
  geom_hline(color = "blue3", yintercept=0) + stat_smooth(se = FALSE, method = "loess", color = "red3")

# Heatmap
mat.dist = pseudoCount
colnames(mat.dist) = paste(colnames(mat.dist), libType, sep = " : ")
mat.dist = as.matrix(dist(t(mat.dist)))
mat.dist = mat.dist/max(mat.dist)
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(16)
cim(mat.dist, color = rev(hmcol), symkey = FALSE, margins = c(9,9))

# PCA
rv = rowVars(as.matrix(pseudoCount))
ntop = 100
select = order(rv, decreasing = TRUE)[1:ntop]
pca = prcomp(t(pseudoCount[select, ]))

annot = AnnotatedDataFrame(data = data.frame(condition, libType,
                                             row.names = colnames(pseudoCount)))
expSet = new("ExpressionSet", exprs = as.matrix(pseudoCount), phenoData = annot)
plotPCA(DESeqTransform(expSet), intgroup = c("condition", "libType"))

pcoa = rda(as.matrix(pseudoCount))
biplot(pcoa)

# MDS
x = pseudoCount
s = rowMeans((x - rowMeans(x))^2)
o = order(s, decreasing = TRUE)
x = x[o, ]
x = x[1:ntop, ]
D[i, j] = sqrt(colMeans((x[, i] - x[, j])^2))

# Filtering
keep = rowSums(count) > 0
filtCount = pseudoCount[keep, ]
dim(count)
