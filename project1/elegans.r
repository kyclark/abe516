source("http://www.bioconductor.org/biocLite.R")
library("affy")
library("simpleaffy")
library("scales")
library("R.utils")

setwd("~/work/abe516/project1")
cwd = getwd()

dat = ReadAffy(celfile.path = "data")
featureNames(dat)[1:10]
num_samples = length(dat)
printf("There are %s features in %s samples\n", len(featureNames(dat)), num_samples)

head(pm(dat, "172682_x_at"))
head(mm(dat, "172682_x_at"))
hist(dat)

snames = c('c1', 'c2', 'bs1', 'bs2', 'bs3', 'efs1', 'efs2', 'efs3', 'efm1', 'efm2', 'efm3', 'c3', 'c4')
groups = paste0("G", unlist(strsplit("3311100022233", '')))
ngroups = len(unique(groups))

groups = c("C", "Bs", "Efs", "Efm")
fl <- factor(c(rep('0', 2), rep('1', 3), rep('2', 3), rep('3', 3), rep('0', 2)), 
             labels=groups)
palette(brewer_pal(type = "seq", palette = "Set2")(len(groups)))
boxplot(dat, names = fl, main = "Raw", las = 2, col = fl)
plot(probeset(dat, geneNames(dat)[1]) [[1]])

RNAdeg = AffyRNAdeg(dat)
plotAffyRNAdeg(RNAdeg)

png(filename = file.path(cwd, "ma-raw.png"), width = 600, height = 800)
par(mfrow = c(5,3))
MAplot(dat, plot.method = "smoothScatter")
dev.off()

dat.rmabg = bg.correct(dat, "rma")
png(filename = file.path(cwd, "ma-bg-corrected.png"), width = 600, height = 800)
par(mfrow = c(5,3))
MAplot(dat.rmabg, plot.method = "smoothScatter")
dev.off()

dat.norm = normalize(dat.rmabg, "quantiles")
png(filename = file.path(cwd, "ma-normalized.png"), width = 600, height = 800)
par(mfrow = c(5,3))
MAplot(dat.norm, plot.method = "smoothScatter")
dev.off()

boxplot(dat.norm, 
        names = snames, 
        main = "Normalized", 
        las = 2,
        col = as.factor(groups))

# Nothing from here on works
# what can I do with this?
expr.dat <- expresso(dat, 
                     bgcorrect.method = "rma", 
                     normalize.method = "constant", 
                     pmcorrect.method = "pmonly", 
                     summary.method   = "avgdiff")

# does not compute
hist(expr.dat)
MAplot(expr.dat, plot.method = "smoothScatter")
#top250 = read.table("top250.txt", header = T)
