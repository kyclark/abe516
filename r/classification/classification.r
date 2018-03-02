# source("http://www.bioconductor.org/biocLite.R")
# biocLite("ALL")
# biocLite("genefilter")
# biocLite("hgu95av2.db")
# biocLite("MLInterfaces")
# install.packages("gplots")
# install.packages("e1071")

library("ALL") 
library("genefilter") 
library("hgu95av2.db") 
library("MLInterfaces") 
library("gplots")
library("e1071")

data("ALL")
exprs(ALL)[1:3,]
pData(ALL)
names(pData(ALL))
ALL.1=ALL[,order(ALL$mol.bio)]
head(ALL.1)
ALL.1$mol.bio

heatmap(cor(exprs(ALL.1)), 
        Rowv = NA, 
        Colv = NA,
        scale = "none", 
        labRow = ALL.1$mol.bio, 
        labCol = ALL.1$mol.bio, 
        RowSideColors = as.character(as.numeric(ALL.1$mol.bio)),
        ColSideColors = as.character(as.numeric(ALL.1$mol.bio)), 
        col = greenred(75))

# Make a simple filtering and select genes with standard deviation 
# (i.e., sd) larger than 1 
ALL.sd = apply(exprs(ALL.1), 1, sd)  

# 379 genes left 
ALL.new = ALL.1[ALL.sd>1, ] 
ALL.new  

heatmap(cor(exprs(ALL.new)), 
        Rowv = NA, 
        Colv = NA,
        scale = "none", 
        labRow = ALL.new$mol.bio, 
        labCol = ALL.new$mol.bio, 
        RowSideColors = as.character(as.numeric(ALL.new$mol.bio)),
        ColSideColors = as.character(as.numeric(ALL.new$mol.bio)), 
        col = greenred(75))

table(ALL.new$BT)

bcell = grep("^B", as.character(ALL$BT))
ALL_B=ALL[, bcell]

table(ALL$mol.biol)
ALL_bcrneg = ALL_B[,ALL_B$mol.biol=="NEG"|ALL_B$mol.biol=="BCR/ABL"]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)

# Apply non-specific filtering
library(genefilter)
small = nsFilter(ALL_bcrneg, var.cutoff=0.75)$eset
small     ## find out how many genes are left?
dim(small) 
table(small$mol.biol)

Negs = which(small$mol.biol == "NEG")
Bcr = which(small$mol.biol == "BCR/ABL")

set.seed(7)
S1 = sample(Negs, 20, replace=FALSE)
S2 = sample(Bcr, 20, replace=FALSE)
TrainInd = c(S1, S2)

Train=small[,TrainInd]

ordTT = order(abs(Traintt$statistic), decreasing=TRUE)
esetShort <- small[fname50,]
dim(esetShort)
