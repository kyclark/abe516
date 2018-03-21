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
ALL_bcrneg = ALL_B[ , ALL_B$mol.biol == "NEG" | ALL_B$mol.biol == "BCR/ABL"]
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

Train = small[, TrainInd]

Traintt = rowttests(Train, "mol.biol") # run two-sample t-test
head(Traintt)

ordTT = order(abs(Traintt$statistic), decreasing = TRUE)
fname50 = featureNames(small[ordTT[1:50], ])
esetShort = small[fname50,]
dim(esetShort)

# K nearest neighbors
Knn.out = MLearn(mol.biol ~ ., data = esetShort, knnI(k=1, l=0), TrainInd)
show(Knn.out)
confuMat(Knn.out)
bb = confuMat(Knn.out)
Err = (bb[1,2] + bb[2,1]) / sum(bb) # calculate misclassification rate
(sum(bb) - sum(diag(bb))) / sum(bb)

# Linear discriminant analysis
Lda.out = MLearn(mol.biol ~ ., data = esetShort, ldaI, TrainInd)
show(Lda.out)
confuMat(Lda.out)
bb = confuMat(Lda.out)
Err = (bb[1,2] + bb[2,1]) / sum(bb) # calculate misclassification rate

(sum(bb) - sum(diag(bb))) / sum(bb)

# Support vector machine
Svm.out = MLearn(mol.biol ~ ., data = esetShort, svmI, TrainInd)
show(Svm.out)
bb = confuMat(Svm.out)
Err = (bb[1,2]+bb[2,1])/sum(bb) # calculate misclassification rate
(sum(bb) - sum(diag(bb))) / sum(bb)

# Classification tree
Ct.out = MLearn(mol.biol ~ ., data = esetShort, rpartI, TrainInd)
show(Ct.out)
bb = confuMat(Ct.out)
(sum(bb) - sum(diag(bb))) / sum(bb)

# Logistic regression
Lr.out = MLearn(mol.biol ~ ., data = esetShort, glmI.logistic(threshold=.5), TrainInd,
                family=binomial)
show(Lr.out)
bb = confuMat(Lr.out)
(sum(bb) - sum(diag(bb))) / sum(bb)


ranpart = function(K, data) {
  N = nrow(data)
  cu = as.numeric(cut(1:N, K))
  sample(cu, size = N, replace = FALSE)
}

ranPartition = function(K) function(data, clab, iternum) {
  p = ranpart(K, data)
  which(p != iternum)
}


set.seed(1)
kk=100
error=rep(0, kk)
for (i in 1:kk){
  r1 = MLearn(mol.biol ~ ., esetShort, knnI(k = 1, l = 0), 
              xvalSpec("LOG", 5, partitionFunc = ranPartition(5)))
  bb= confuMat(r1)
  error[i]=(bb[1,2]+bb[2,1])/sum(bb)
}
mean(error) ## get the average of the error/ misspecification from cross-validation
boxplot(error)
