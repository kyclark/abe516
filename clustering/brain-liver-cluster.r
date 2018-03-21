library(affy)
library(limma) ##– linear model for DE
library(Heatplus)## – for heat map view
library(gplots) ## - for color change in heatmap plot

data.dir = "/Users/kyclark/work/abe516/r/de/data"
cel.files = file.path(data.dir, dir(data.dir, pattern = "*.CEL"))
raw = ReadAffy(filenames = cel.files)
eset = rma(raw)

# construct a design matrix for this experiment
f = factor(c(1, 1, 2, 2, 3, 3, 4, 4), 
           labels=c("brain", "f.brain", "f.liver", "liver"))
design = model.matrix(~ 0 + f)
colnames(design) = c("brain", "f.brain", "f.liver", "liver")
fit = lmFit(eset, design)
contrast.matrix = makeContrasts(f.brain-brain, f.liver - liver, 
                                levels = design)
fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)

## select genes with adjusted p-value < 0.0001
p.bh = p.adjust(fit2$F.p.value, method = "BH")
select = (p.bh < 0.0001)
table(select) ## how many genes are selected for further analysis

newdata = exprs(eset)[select,]
#heatmap_2(newdata, legend = 4)
heatmap_2(newdata, legend = 4, col = greenred(20), legfrac = 10)

# Change the grouping of columns:
# create a hierarchical tree "by hand" and cut it
hc = hclust(dist(t(newdata)))
c4 = cutree(hc, k = 4)
heatmap_plus(newdata, clus = c4, addvar = design, col = greenred(20))

c2 = cutree(hc, k = 2)
heatmap_plus(newdata, addvar = design, clus = c2, col = greenred(20))

# "Scale" parameter of heatmap
# default value is “row”
# only impacts the colors of the image. # It does not affect the dendrograms.
heatmap_2(newdata, scale = "row", legend = 4, col = greenred(20)) ### we already saw it
heatmap_2(newdata, scale = "none", legend = 4, col = greenred(20))

# But scaling the data ...
# These two codes generate different dendrograms

# If you focus on the absolute distance among genes, use:
heatmap_2(newdata, scale="row", legend=4,col=greenred(20))  ### we already saw it

# If you’re interested in the expression pattern of genes across experiments, use
heatmap_2(t(scale(t(newdata))), scale="none", legend=4, col=greenred(20) )
