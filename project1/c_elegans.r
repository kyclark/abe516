#!/usr/bin/env Rscript

library("R.utils")
source("http://www.bioconductor.org/biocLite.R")
library("affy")
library("simpleaffy")
library("scales")
library("limma")
library("SAGx")
library("Hmisc")
library("Heatplus")
library("gplots")
library("cluster")
library("GEOquery")
library("topGO")
library("gridExtra")

#cwd = getwd() # this works from the command line but not in RStudio
cwd = "/Users/kyclark/work/abe516-project1"
setwd(cwd)
data_dir = file.path(cwd, 'data')
if (!dir.exists(data_dir)) {
  stop(paste0("Missing data dir ", data_dir))
}

figures_dir = file.path(cwd, 'figures')
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir)
}

tables_dir = file.path(cwd, 'tables')
if (!dir.exists(tables_dir)) {
  dir.create(tables_dir)
}

dat = ReadAffy(celfile.path = "data")
num_samples = length(dat)
groups = c("C", "Bs", "Efs", "Efm")
fl <- factor(c(rep('0', 2), rep('1', 3), rep('2', 3), rep('3', 3), rep('0', 2)), 
             labels=groups)

printf("There are %s features in %s samples\n", length(featureNames(dat)), num_samples)

#
# Raw plots
# 
printf("Examing raw values\n")
palette(brewer_pal(type = "seq", palette = "Set2")(length(groups)))
png(filename = file.path(figures_dir, "raw-box-plot.png"))
boxplot(dat, names = fl, main = "Raw", las = 2, col = fl)
invisible(dev.off())

png(filename = file.path(figures_dir, "raw-ma-plot.png"), width = 600, height = 800)
par(mfrow = c(5,3))
MAplot(dat, plot.method = "smoothScatter")
invisible(dev.off())

dat.rmabg = bg.correct(dat, "rma")
dat.norm = normalize(dat.rmabg, "quantiles")
png(filename = file.path(figures_dir, "norm-box-plot.png"))
boxplot(dat.norm, main = "Normalized", las = 2, col = as.factor(groups))
invisible(dev.off())

png(filename = file.path(figures_dir, "norm-ma-plot.png"), width = 600, height = 800)
par(mfrow = c(5,3))
MAplot(dat.norm, plot.method = "smoothScatter")

#
# Rather than manually running each error-correcting step, we can use "expresso"
#
print("Running expresso")
eset <- expresso(dat, 
                 bgcorrect.method = "rma", 
                 normalize.method = "quantiles", 
                 pmcorrect.method = "pmonly", 
                 summary.method   = "medianpolish")

design <- model.matrix(~ 0 + fl) # original
colnames(design) = groups
#
# To get the metadata we want to see in the output, it's necessary to use the SOFT file
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95636/soft/GSE95636_family.soft.gz
# The GEOquery library provides a "getGEO" function that will download this and 
# provide us with the corrected expression data.
#
gset <- getGEO("GSE95636", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep("GPL200", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)

fit <- lmFit(gset, design)
contrast.matrix <- makeContrasts(Efm-C, Efs-C, Bs-C, Efs-Bs, Efm-Efs, levels = design)
fit1 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit1)
tt.all = topTable(fit2, 
                  adjust = "fdr", 
                  p.value = 0.05, 
                  lfc = 1, 
                  sort.by = "B", 
                  number = 250)
write.table(tt.all, 
            file = file.path(tables_dir, "top_overall.tab"), 
            row.names = FALSE, 
            sep = "\t")

for (coef in 1:5) {
  filename = file.path(tables_dir, paste0("top250-c", coef, ".tab"))
  printf("Writing '%s'\n", filename)
  tt = topTable(fit2, 
                coef = coef, 
                adjust = "fdr", 
                sort.by = "B", 
                lfc = 1,
                number = 250)
  
  write.table(tt, file = filename, row.names = F, sep = "\t")  
}

results = decideTests(fit2, 
                      adjust = "fdr", 
                      p = 0.05, 
                      lfc = 1)

summary(results)

printf("Writing Venn diagrams\n")
png(filename = file.path(figures_dir, "venn.png"), width = 800, height = 400)
par(mfrow = c(1, 2))
#vennDiagram(results, main = "Diff Exp Genes")
vennDiagram(results, include = "down", main = "Down Regulated")
vennDiagram(results, include = "up", main = "Up Regulated")
invisible(dev.off())

printf("Writing volcano plots\n")
png(filename = file.path(figures_dir, "volcano.png"), 
    width = 800, 
    height = 800)
par(mfrow = c(2,3))
volcanoplot(fit2, 
            coef = 1, 
            main = "E_faecium-Control", 
            names = row.names(fit2$coefficients), 
            highlight = 10)
volcanoplot(fit2, 
            coef = 2, 
            main = "B_subtilis-Control", 
            names = row.names(fit2$coefficients), 
            highlight = 10)
volcanoplot(fit2, 
            coef = 3, 
            main = "E_faecium-B_subtilis", 
            names = row.names(fit2$coefficients), 
            highlight = 10)
volcanoplot(fit2, 
            coef = 4, 
            main = "E_faecium-E_faecalis", 
            names = row.names(fit2$coefficients), 
            highlight = 10)
volcanoplot(fit2, 
            coef = 5, 
            main = "E_faecium-E_faecalis", 
            names = row.names(fit2$coefficients), 
            highlight = 10)
invisible(dev.off())

#
# Find the top GO functions; many of the function fields have multiple
# functions separated by "///", so we split on that; additionally, some
# of the function names are very long, so we remove (gsub) everything after
# a comma to get a usable name. As there are so many, we use only those 
# occuring with a frequency greater than X
#
go_functions = sub(",.*", "", 
                   unlist(strsplit(tt.all$GO.Function, '///', fixed = T)))
function_count = as.data.frame(table(go_functions))
colnames(function_count) = c("func", "freq")
top_functions = function_count[function_count$freq > 3,]

#
# There are too many processes to show, so we need to filter 
# to those with a frequency greater than X
#
go_processes = unlist(strsplit(tt.all$GO.Process, '///', fixed = T))
process_count = as.data.frame(table(go_processes))
colnames(process_count) = c("process", "freq")
top_processes = process_count[process_count$freq > 3,]

#
# This code is ugly, but it's necessary to force ggplot to keep the order
#
top_functions$func = factor(top_functions$func, 
                            levels = top_functions$func[order(top_functions$freq)])
top_processes$process = factor(top_processes$process, 
                               levels = top_processes$process[order(top_processes$freq)])

#
# ggplot doesn't play nicely with par, so we have to use gridExtra
#
png(filename = file.path(figures_dir, "go.png"), 
    width = 800,
    height = 600)
p1 = ggplot(data=top_functions, aes(x=func, y=freq)) + 
  geom_bar(stat="identity", fill="#56B4E9", colour="black") + 
  coord_flip()
p2 = ggplot(data=top_processes, aes(x=process, y=freq)) + 
  geom_bar(stat="identity", fill="#56B4E9", colour="black") + 
  coord_flip()
grid.arrange(p1, p2, ncol=2, nrow = 1)
invisible(dev.off())

#
# Genomic location
#
locs = unlist(strsplit(tt.all$Chromosome.annotation, '///', fixed = T))
chromosome = sub('Chromosome ', '', sub(',.*', '', locs))
png(filename = file.path(figures_dir, 'chromosomes.png'))
qplot(chromosome)
invisible(dev.off())

#
# Clustering
# 
printf("Cluster analysis\n")
#get the expression matrix form eset
eset.expr = exprs(eset)
#extract overall differential expressed genes from expression matrix
DE.all = eset.expr[rownames(tt.all),]
#compute distance matrix
dist.cor = as.dist(1 - cor(t(DE.all)))
#hierarchical clustering
hcfit = hclust(dist.cor, method = 'ave')

#
# Using gap statistic to determine k in HC
#
png(filename = file.path(figures_dir, "gapstat.png"), 
    width = 800, 
    height = 800)
k = 6
Gap = rep(0,k)
se = rep(0,k)
for (i in 2:k) {
  mem = cutree(hcfit, i)
  result = gap(DE.all, class=mem)
  Gap[i] = result[1]
  se[i] = result[2]
}
errbar(1:k, Gap, Gap-se, Gap+se, xlab = "Number of clusters")
lines(1:k, Gap) 
invisible(dev.off())

#
# Using silhoutte width to determine k
#
png(filename = file.path(figures_dir, "silhouette.png"), 
    width = 800, 
    height = 800)
k = 10
sil = rep(0,k)
for (i in 2:k){
  mem = cutree(hcfit,i)
  aa = silhouette(mem, dist(DE.all))
  sil[i] = mean(aa[,3])
}
plot(1:k,sil)
lines(1:k,sil)
invisible(dev.off())

#
# Heatmap
#
png(filename = file.path(figures_dir, "heatmap.png"), 
    width = 800,
    height = 800)
c2 = cutree(hcfit,k=2) 
colnames(DE.all) = c(rep('E. coli', 2), 
                     rep('B. subtilis', 3), 
                     rep('E. faecalis', 3), 
                     rep('E. faecium', 3),
                     rep('E. coli', 2))
heatmap_plus(t(scale(t(DE.all))), 
             clus = c2, 
             addvar = design, 
             col = greenred(20))
invisible(dev.off())

#
# plot the resulting dendrogram
#
png(filename = file.path(figures_dir, "dendrogram.png"), 
    width = 800, 
    height = 400)
plot(hcfit)
invisible(dev.off())

# cut the tree at height 1.25
hcfit1 = cutree(hcfit, h = 1.25)

# generate 2 clusters and save them in a table
clus = table(names(hcfit1), hcfit1)
write.table(clus, 
            file = file.path(tables_dir, "clus.tab"), 
            row.names = TRUE, 
            sep = "\t")
invisible(dev.off())

printf("All done!\n")