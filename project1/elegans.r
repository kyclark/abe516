source("http://www.bioconductor.org/biocLite.R")
library("affy")
library("simpleaffy")
library("scales")
setwd("~/work/abe516/project1")
cwd = getwd()

dat = ReadAffy(celfile.path = "data")
featureNames(dat)[1:10]
num_samples = length(dat)

head(pm(dat, "172682_x_at"))
head(mm(dat, "172682_x_at"))
hist(dat)

snames = c('c1', 'c2', 'bs1', 'bs2', 'bs3', 'efs1', 'efs2', 'efs3', 'efm1', 'efm2', 'efm3', 'c3', 'c4')
groups = paste0("G", unlist(strsplit("3311100022233", '')))
ngroups = len(unique(groups))
palette(brewer_pal(type = "seq", palette = "Set2")(ngroups))
boxplot(dat, names = snames, main = "Raw", las = 2, col = as.factor(groups))
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
