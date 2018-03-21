library(SAGx) ## for gap statistic
library(cluster) ## for hierarchical cluster, k-mean cluster & silhouette width
library(Hmisc) ## for error bar plot

set.seed(7)       ## for random number later
x1 <- c(rnorm(20, sd = .05), 
        rnorm(20, mean = 1, sd = .05),
        rnorm(20, mean = 2.5, sd= .05))

## x1 dimension with 3 clusters
x2 <- 4+c(rnorm(20, sd = 0.5),
          rnorm(20, mean = 10, sd = 1.0), 
          rnorm(20, mean = 10, sd = 1))

## x2 dimension with 3 clusters
plot(x1,x2)  ## scatter plot of the data

dat = cbind(x1, x2) ##  combine two vectors into one data matrix
hc0 = hclust(dist(dat), method = "ave")
## hierarchical clustering (hc) with average linkage
plot(hc0) ## plot the dendrogram of the hc result

# Using gap statistic to determine k in HC - 1
k = 6   ## check all possible numbers of clusters, 2~6
gap = rep(0, k) ## initialize the gap statistics
se = rep(0, k) ## initialize the standard error
for (i in 2:k) {
  mem = cutree(hc0, i) ## get the cluster membership by using “cuttree” on the object of hier. cluster. 
  result=gap(dat, class=mem) 
  ## get the gapstatistics
  gap[i] = result[1] ## extract the gap stat values se[i]=result[2] ## and the s.e. values
}

errbar(1:k, gap, gap - se, gap + se, xlab = "Number of clusters") ## error bar plot
lines(1:k, gap) ## connect them

# Using silhouette width to determine k in HC
k = 6
sil = rep(0,k)
for (i in 2:k){
  mem = cutree(hc0, i)
  aa = silhouette(mem, dist(dat))
  sil[i] = mean(aa[,3])
}
plot(1:k, sil)
lines(1:k, sil)

km = kmeans(dat, centers = 2)
mem = km$cluster
plot(x1, x2, col = mem+1)

# Use gap statistic to determine k in K-means
k = 6   ## check all possible numbers of clusters, 2~6
gap = rep(0,k)      ## initialize the gap statistics
se = rep(0, k) ## initialize the standard error
for (i in 2:k) {
  km = kmeans(dat, centers = i)
  mem = km$cluster
  result = gap(dat, class = mem) ## get the gap statistics
  gap[i] = result[1] ## extract the gap stat values se[i]=result[2] ## and the s.e. values
}

errbar(1:k, gap, gap-se, gap+se, xlab="Number of clusters") ## error bar plot
lines(1:k, gap) ## connect them

# Use silhouette width to determine k in K-means
k = 6
sil = rep(0, k)
for (i in 2:k) {
  km = kmeans(dat, centers = i)
  mem = km$cluster
  aa = silhouette(mem, dist(dat))
  sil[i] = mean(aa[,3])
}
plot(1:k, sil)
lines(1:k, sil)

# How about just let k=3 for both HC and K- means, and check
hc0 = hclust(dist(dat))
mem = cutree(hc0, 3)
plot(x1, x2, col = mem)

km = kmeans(dat, 3)
mem = km$cluster
plot(x1, x2, col = mem)

# revisit the data ...
plot(x1, x2)
plot(x1, x2, xlim = c(-0.2, 16), ylim = c(-0.2, 16))

# Scaling needed
x1.sc = x1/sd(x1)
x2.sc = x2/sd(x2)
dat.sc = cbind(x1.sc, x2.sc)
hc1 = hclust(dist(dat.sc), "ave")
plot(hc1)

# For the scaled data, perform HC and gap
K = 6
gap = rep(0, k)
se = rep(0, k)

for (i in 2:k) {
  mem = cutree(hc1, i)
  result = gap(dat.sc, class = mem)
  gap[i] = result[1]
  se[i]=result[2]
}
errbar(1:k, gap, gap-se, gap+se, xlab="Number of clusters")
mem3 = cutree(hc1, 3)
plot(x1, x2, col = mem3)
lines(1:k, gap)

# For the scaled data, perform Kmeans and gap 
k = 6
gap = rep(0,k)
se = rep(0,k)
for (i in 2:k) {
  km = kmeans(dat.sc, centers = i)
  mem = km$cluster
  result = gap(dat.sc, class = mem)
  gap[i] = result[1]
  se[i] = result[2]
}

errbar(1:k, gap, gap-se, gap+se, xlab="Number of clusters")
lines(1:k, gap)

# For the scaled data, perform Kmeans and silhouette width
k = 6
sil = rep(0, k)
for (i in 2:k) {
  km = kmeans(dat.sc, centers=i)
  mem = km$cluster
  aa = silhouette(mem, dist(dat.sc))
  sil[i] = mean(aa[,3])
}
plot(1:k, sil)
lines(1:k, sil)

km = kmeans(dat.sc, centers = 3)
mem3 = km$cluster
# What???? K-means is sensitive to seeds initialization
plot(x1, x2, col=mem3)

# Recommend: Hybrid ... (hierarchical+kmeans)
# Use the seeds info from HC result
hc = hclust(dist(dat.sc), method = "ave")
mem = cutree(hc, 3)
c1 = tapply(x1.sc, mem, mean)
c2 = tapply(x2.sc, mem, mean)

# then do kmeans clustering:
km = kmeans(dat.sc, centers = cbind(c1, c2))
plot(x1, x2, col=km$cluster)
