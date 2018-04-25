library("ggplot2")
library("ggdendro")
library("pheatmap")
library("RColorBrewer")
library("R.utils")

wd = "~/work/abe516/project3"
setwd(wd)
full.df = read.table("dist.txt", header = T, check.names = F)
colors = colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(full.df, col = colors, filename = file.path(wd, "full-heatmap.png"))

# Subsample full dataset
num.samples = 50
num.iter = 10
for (i in 1:num.iter) {
  print(paste("Sample", i))
  sample.cols = sample(nrow(full.df), num.samples)
  df = full.df[sample.cols, sample.cols]
  fit = hclust(as.dist(df), method = "ward.D2")
  dg = ggdendro::ggdendrogram(fit, rotate=F)
  ggsave(file = file.path(wd, paste0("sample", i, ".png")),
         limitsize = FALSE, width = 10, height = 5, plot = dg)
}

# HOT
hot.df = read.table("hot-dist.txt", header = T, check.names = F)
hot.fit = hclust(as.dist(hot.df), method = "ward.D2")
dg = ggdendro::ggdendrogram(hot.fit, rotate=F)
ggsave(file = file.path(wd, "hot.png"),
       limitsize = FALSE, width = 10, height = 5, plot = dg)
pheatmap(hot.df, 
         col = colors, 
         cutree_rows = 2,
         cutree_cols = 2,
         filename = file.path(wd, "hot-heatmap.png"))