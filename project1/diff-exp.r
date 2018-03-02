source("http://www.bioconductor.org/biocLite.R")
#biocLite("affy")
#biocLite("limma")
#biocLite("hgu95acdf")

library("affy")
library("limma")
#library("hgu95acdf") 

setwd("~/work/abe516/project1")
cwd = getwd()
raw = ReadAffy(celfile.path = "data")
probeNames(raw)[1:20]
geneNames(raw)[1:10]
#image(raw[,1])
hist(raw)
#boxplot(raw, col=c(2,2,3,3,4,4,5,5))
#par(mfrow=c(2,4))
#MAplot(raw, plot.method="smoothScatter")

# raw.bg=bg.correct(raw, "rma")
# raw.norm=normalize(raw.bg, "quantiles")
# hist(raw.norm, main = " raw data after RMA  BG & quantile normalization ")
# par(mfrow=c(2,4))
# MAplot(raw.norm, plot.method= "smoothScatter")
# raw.pm=pmcorrect.pmonly(raw.norm)
# eset <- computeExprSet(raw.norm, "pmonly","medianpolish")

eset <- expresso(raw, 
                 bgcorrect.method = "rma", 
                 normalize.method = "quantiles", 
                 pmcorrect.method = "pmonly", 
                 summary.method   = "medianpolish")
groups = c("C", "Bs", "Efs", "Efm")
fl <- factor(c(rep('0', 2), rep('1', 3), rep('2', 3), rep('3', 3), rep('0', 2)), 
             labels=groups)
design <- model.matrix(~ 0 + fl)
colnames(design) = groups

fit <- lmFit(eset, design)

contrast.matrix <- makeContrasts(Efm-C, Bs-C, Efs-Bs, Efm-Efs, levels = design)

fit1 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit1)

topTable(fit2, n = 10, adjust = "fdr")
topTable(fit2, coef = 1, n = 10, adjust = "fdr")
topTable(fit2, coef = 2, n = 10, adjust = "fdr")

modt = topTable(fit2, 
                coef    = 1,
                sort.by = "t",
                adjust  = "fdr",
                number  = 10)[,c(1,3)]

logodds = topTable(fit2, 
                   coef    = 1,
                   sort.by = "B",
                   adjust  = "fdr",
                   number  = 10)[,c(1,6)]

logFC = topTable(fit2, 
                 coef    =1,
                 sort.by = "logFC",
                 adjust  = "fdr",
                 number  = 10)[,c(1,3)]

results  = decideTests(fit2, adjust = "fdr", p = 0.05)
results2 = decideTests(fit2, adjust = "fdr", p = 0.005)

summary(results)
summary(results2)
p.value <-fit2$F.p.value

# i <- grep("AFFX", geneNames(raw))
# summary(p.value[i])

results <- classifyTestsF(fit2, p.value = 0.000001)
head(results)
summary(results)
vennDiagram(results)

table(efm.c=results[,1], bs.c=results[,2])

png(filename = file.path(cwd, "venn.png"), width = 400, height = 800)
par(mfrow = c(3, 1))
vennDiagram(results, main = "Diff Exp Genes")
vennDiagram(results, main = "Down", include = "down")
vennDiagram(results, main = "Up", include = "up")
dev.off()

png(filename = file.path(cwd, "volcano2.png"), width = 800, height = 800)
par(mfrow = c(2,2))
volcanoplot(fit2,
            coef=1, 
            main="E_faecium-Control",
            names=row.names(fit2$coefficients), 
            highlight=10)

volcanoplot(fit2,
            coef=2, 
            main="B_subtilis-Control",
            names=row.names(fit2$coefficients), 
            highlight=10)

volcanoplot(fit2,
            coef=3, 
            main="E_faecium-B_subtilis",
            names=row.names(fit2$coefficients), 
            highlight=10)

volcanoplot(fit2,
            coef=3, 
            main="E_faecium-E_faecalis",
            names=row.names(fit2$coefficients), 
            highlight=10)
dev.off()

# brain.up <- results[results[,1]==1, 1]
# brain.down <- results[results[,1]== -1,1 ]
# write.csv(data.frame(brain.up),file="myresults_brain_up.csv")
# write.csv(data.frame(brain.down),file="myresults_brain_down.csv")
# 
# par(mfrow=c(1,2))
# volcanoplot(fit2,coef=1, main="brain",
#             names=row.names(fit2$coefficients), highlight=10)
# volcanoplot(fit2,coef=2, main="liver",
#             names=row.names(fit2$coefficients), highlight=10)
