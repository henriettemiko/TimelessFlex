
##########
#name:          plot_clusters_coherence_fold.r
#description:   plots cluster coherence for fold changes subset
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          November 21, 2019
##########


library(psych)

setwd(getwd())

args = commandArgs(trailingOnly=TRUE)

numClusters = as.numeric(args[1])
numMarks = as.numeric(args[2])
numTimePoints = as.numeric(args[3])
numcluster.dir = args[4]
timeless.dir = args[5]
signalgenerator.dir = args[6]
type = args[7]

score <- function(x){
    x = scale(x, scale = FALSE)
    x  = ((x-min(x))/(max(x)-min(x))) * 100
    return(x)
}

l = read.table(paste0(numcluster.dir, "/classes-", numClusters, ".txt"))
k = read.table(paste0(signalgenerator.dir, "/allFold_data.txt"))

# 3 fold changes
k27ac = cbind(k[[1]], k[[2]], k[[3]])
k27me3 = cbind(k[[4]], k[[5]], k[[6]])
k4me1 = cbind(k[[7]], k[[8]], k[[9]])
k4me3 = cbind(k[[10]], k[[11]], k[[12]])

k27ac = matrix(score(as.vector(k27ac)), ncol = 3, byrow = FALSE)
k27me3 = matrix(score(as.vector(k27me3)), ncol = 3, byrow = FALSE)
k4me1 = matrix(score(as.vector(k4me1)), ncol = 3, byrow = FALSE)
k4me3 = matrix(score(as.vector(k4me3)), ncol = 3, byrow = FALSE)

print(colMeans(k27ac))
print(colMeans(k27me3))
print(colMeans(k4me1))
print(colMeans(k4me3))

cbbPalette <- c("gray", "#009E73", "#D50F25", "black")
#me3 gray, me1 black, 27ac green 009E73, 27me3 red
#plot from dark to light colors: me1, k27me3, k27ac, me3


coherence <- rep(-1,numClusters)
correlation.means <- rep(-1,numClusters)

for (i in 1:numClusters) {


    cur.items = which(l[[1]] == i)
    stuff.k27ac <- cbind(k27ac[cur.items,1], k27ac[cur.items,2], 
                         k27ac[cur.items,3]) 
    stuff.k27me3 <- cbind(k27me3[cur.items,1], k27me3[cur.items,2], 
                          k27me3[cur.items,3]) 
    stuff.k4me1 <- cbind(k4me1[cur.items,1], k4me1[cur.items,2], 
                         k4me1[cur.items,3]) 
    stuff.k4me3 <- cbind(k4me3[cur.items,1], k4me3[cur.items,2], 
                         k4me3[cur.items,3])

    stuff.all <- cbind(stuff.k27ac,stuff.k27me3,stuff.k4me1,stuff.k4me3)

    mean.all <- colMeans(stuff.all)

    cor.all <- apply(stuff.all,1,cor,mean.all)
    #cor.all <- apply(stuff.all,1,function(r) cor(r,mean.all, method="spearman"))

    cor.all.len <- length(cor.all)

    #all with correlation bigger than 0.8
    cor.all.good <- cor.all[cor.all>=0.8]
    cor.all.good.len <- length(cor.all.good)

    #cluster coherece:
    #fraction of regions with a Pearson correlation of at least 0.8 to
    #the cluster mean (like CMINT)
    coherence[i] = cor.all.good.len/cor.all.len

    str(cor.all)

    print(summary(cor.all))

    correlation.means[i] <- mean(cor.all)
}


print(coherence)
print(correlation.means)

write.table(coherence,file=paste0(type, "_", numClusters, "_coherence_fold.txt"), sep="\t", quote=F, row.names=F, col.names=F)

write.table(correlation.means,file=paste0(type, "_", numClusters, "_correlations_fold.txt"), sep="\t", quote=F, row.names=F, col.names=F)

q()


pdf(paste0(type, "_clusters_", numClusters, "_coherence_fold.pdf"), height=6, width=8)
plot(coherence, xlab="clusters", ylab="cluster coherence", ylim=c(0,1), pch=19)
axis(side=1, at=1:numClusters, labels=1:numClusters)
dev.off()

pdf(paste0(type, "_clusters_", numClusters, "_coherence_correlations_fold.pdf"), height=6, width=8)
plot(correlation.means, xlab="clusters", ylab="cluster coherence", ylim=c(0,1), pch=19)
axis(side=1, at=1:numClusters, labels=1:numClusters)

dev.off()
warnings()


q()

