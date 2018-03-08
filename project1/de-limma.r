# Differential expression analysis with limma
source("http://www.bioconductor.org/biocLite.R")
library("Biobase")
library("GEOquery")
library("limma")

# load series and platform data from GEO
setwd("~/work/abe516/project1")
cwd = getwd()

gset <- getGEO("GSE95636", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep("GPL200", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "3311100022233"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
groups = c("C", "Bs", "Efs", "Efm")
fl <- factor(c(rep('0', 2), rep('1', 3), rep('2', 3), rep('3', 3), rep('0', 2)), 
             labels=groups)
# sml <- paste("G", sml, sep="")    # set group names
# fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
#cont.matrix <- makeContrasts(G3-G0, G1-G0, G2-G1, G3-G2, levels=design)
cont.matrix <- makeContrasts(Efm-C, Bs-C, Efs-Bs, Efm-Efs, levels=design)
#cont.matrix <- makeContrasts(E_faecium-Control, B_subtilis-Control, E_faecalis-B_subtilis, E_faecium-E_faecalis, levels=design)
fit2   = contrasts.fit(fit, cont.matrix)
fit2   = eBayes(fit2, 0.01)
tT     = topTable(fit2, adjust = "fdr", sort.by = "B", number = 250)
cnames = c("ID","adj.P.Val","P.Value","Gene.symbol","Gene.title","GO.Function","GO.Process","GO.Component")
tT      = subset(tT, select = cnames)

write.table(tT, file=file.path(cwd, "top250.tab"), row.names=F, sep="\t")

results = decideTests(fit2, adjust="fdr", p=0.05)
summary(results)
table(efm.c=results[,1], bs.c=results[,2])

png(filename = file.path(cwd, "venn.png"), width = 400, height = 800)
par(mfrow = c(3, 1))
vennDiagram(results, main="Diff Exp Genes")
vennDiagram(results, include="down", main="Down")
vennDiagram(results, include="up", main="Up")
dev.off()

png(filename = file.path(cwd, "volcano.png"), width = 800, height = 800)
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