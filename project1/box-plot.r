# Boxplot for selected GEO samples
source("http://www.bioconductor.org/biocLite.R")
library("Biobase")
library("GEOquery")
library("scales")

# load series and platform data from GEO
setwd("~/work/abe516/project1")
gset <- getGEO("GSE95636", GSEMatrix = TRUE, getGPL = FALSE, destdir = getwd())
if (length(gset) > 1) idx <- grep("GPL200", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
groups     = paste0("G", unlist(strsplit("3311100022233", '')))
ex         = exprs(gset)[, order(groups)] # order samples by group
fl         = as.factor(sort(groups))
labels     = c("E_faecalis", "B_subtilis", "E_faecium", "Control")
num_groups = len(labels)

palette(brewer_pal(type = "seq", palette = "Set2")(num_groups))
dev.new(width = 4 + dim(gset)[[2]]/5, height = 6)
par(mar = c(2 + round(max(nchar(sampleNames(gset))) / 2), 4, 2, 1))
title <- paste ("GSE95636", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, 
        boxwex  = 0.6, 
        notch   = TRUE, 
        main    = title, 
        outline = FALSE, 
        las     = 2, 
        col     = fl)
legend("topleft", labels, fill = palette(), bty="n")