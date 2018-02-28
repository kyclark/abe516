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

