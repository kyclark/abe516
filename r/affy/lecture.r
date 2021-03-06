#source("http://www.bioconductor.org/biocLite.R")
#biocLite("Biobase")
#biocLite("affy")
library("affy")
setwd("~/work/abe516/r/affy/data/")

harvard.rawData = ReadAffy()
geneNames(harvard.rawData)[1:10]
pm(harvard.rawData, "100_g_at")
mm(harvard.rawData, "100_g_at")
hist(harvard.rawData, main = "Harvard data")
boxplot(harvard.rawData, col=c(2,2,2,3,3,3))
plot(probeset(harvard.rawData, geneNames(harvard.rawData)[1]) [[1]])
RNAdeg<-AffyRNAdeg(harvard.rawData)
plotAffyRNAdeg(RNAdeg)

# biocLite("simpleaffy")
library(simpleaffy)
h.qc=qc(harvard.rawData)
avbg(h.qc)
percent.present(h.qc)

MAplot(harvard.rawData, plot.method="smoothScatter", pair=TRUE)
par(mfrow=c(2,3))
MAplot(harvard.rawData, plot.method="smoothScatter")

harvard.exprData <- expresso(harvard.rawData, 
                             bgcorrect.method = "rma", 
                             normalize.method = "constant", 
                             pmcorrect.method = "pmonly", 
                             summary.method = "avgdiff")
