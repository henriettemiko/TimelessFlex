##########
#name:          baseline_kmeans_Zscore.sh
#description:   basline k-means with z-score of unnormalized counts
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          November 21, 2019
##########


library(factoextra)
library(NbClust)
library(psych)
require("cluster")


setwd(getwd())

args = commandArgs(trailingOnly=TRUE)

numClusters = as.numeric(args[1])
numMarks = as.numeric(args[2])
numTimePoints = as.numeric(args[3])
numcluster.dir = args[4]
timeless.dir = args[5]
signalgenerator.dir = args[6]
type = args[7]
name = args[8]


# all counts (before normalization)
#when using this there are no corresponding norm regions for plotting (because of filtering)
#plots are made with these counts directly

data1 = read.table(paste0(signalgenerator.dir, "/all_max_counts_", name, ".txt"))

data1.counts = data1[,-(1:11)]
#z-score transformation
data1.scaled <- scale(data1.counts, center=TRUE, scale=TRUE)


data1.kmeans <- kmeans(data1.scaled, numClusters, iter.max=50, nstart = 5)

fviz_cluster(data1.kmeans, data1.counts, ellipse.type = "norm")+
    theme_minimal()
ggsave(paste0(type, "_", numClusters, "_baseline_Zscore_clus.pdf"), height=6, width=8)

sil <- silhouette(data1.kmeans$cluster, dist(data1.scaled))
fviz_silhouette(sil)

ggsave(paste0(type, "_", numClusters, "_baseline_Zscore_sil.pdf"), height=6, width=8)

write.table(data1.kmeans$cluster, sep="\t", file=paste0(type, "_", numClusters, "_baseline_Zscore.txt"), quote=F, row.names=F, col.names=F)

#######################
#plot

score <- function(x){
    x = scale(x, scale = FALSE)
    x  = ((x-min(x))/(max(x)-min(x))) * 100
    return(x)
}


l = read.table(paste0(numcluster.dir, "/", numClusters, "/", type, "_", numClusters, "_baseline_Zscore.txt"))

# 5 time points
k27ac = data1.counts[,1:5]
k27me3 = data1.counts[,6:10]
k4me1 = data1.counts[,11:15]
k4me3 = data1.counts[,16:20]

k27ac = matrix(score(as.vector(k27ac)), ncol = 5, byrow = FALSE)
k27me3 = matrix(score(as.vector(k27me3)), ncol = 5, byrow = FALSE)
k4me1 = matrix(score(as.vector(k4me1)), ncol = 5, byrow = FALSE)
k4me3 = matrix(score(as.vector(k4me3)), ncol = 5, byrow = FALSE)

cbbPalette <- c("gray", "#009E73", "#D50F25", "black")
#me3 gray, me1 black, 27ac green 009E73, 27me3 red
#plot from dark to light colors: me1, k27me3, k27ac, me3

pdf(paste0(type, "_clusters_", numClusters, "_baseline_Zscore.pdf"), height=6, width=4)

for (i in 1:numClusters) {

    cur.items = which(l[[1]] == i)

    stuff.k27ac <- cbind(k27ac[cur.items,1], k27ac[cur.items,2],
                         k27ac[cur.items,3], k27ac[cur.items,4],
                         k27ac[cur.items,5])
    stuff.k27me3 <- cbind(k27me3[cur.items,1], k27me3[cur.items,2],
                          k27me3[cur.items,3], k27me3[cur.items,4],
                          k27me3[cur.items,5])
    stuff.k4me1 <- cbind(k4me1[cur.items,1], k4me1[cur.items,2],
                         k4me1[cur.items,3], k4me1[cur.items,4],
                         k4me1[cur.items,5])
    stuff.k4me3 <- cbind(k4me3[cur.items,1], k4me3[cur.items,2],
                         k4me3[cur.items,3], k4me3[cur.items,4],
                         k4me3[cur.items,5])

    colnames(stuff.k4me1) <- c("D0", "D2", "D5", "D7", "D10")
    error.bars(stuff.k4me1, eyes = FALSE, sd = FALSE, bars = FALSE,
               arrow.col = cbbPalette[4], labels = colnames(stuff.k4me1),
               ylim = c(0,65),
               ylab = "average histone mark signal",
               xlab="time", col = cbbPalette[4],
               main = paste0("Cluster ", i, ": n=", length(cur.items)))


    lines(colMeans((cbind(k4me1[cur.items,1], k4me1[cur.items,2],
                          k4me1[cur.items,3], k4me1[cur.items,4],
                          k4me1[cur.items,5]))), col = cbbPalette[4],
          ylim = c(0,65), lwd = 8, type="o", pch=16)

    error.bars(stuff.k4me3, eyes = FALSE, sd = FALSE, bars = FALSE,
               arrow.col = cbbPalette[1], ylim = c(0,65),
               col = cbbPalette[1], add = TRUE)
    lines(colMeans((cbind(k4me3[cur.items,1], k4me3[cur.items,2],
                          k4me3[cur.items,3], k4me3[cur.items,4],
                          k4me3[cur.items,5]))), col = cbbPalette[1],
          ylim = c(0,65), lwd = 8, type="o", pch=16)

    error.bars(stuff.k27ac, eyes = FALSE, sd = FALSE, bars = FALSE,
               arrow.col = cbbPalette[2], ylim = c(0,65),
               col = cbbPalette[2], add = TRUE)
    lines(colMeans((cbind(k27ac[cur.items,1], k27ac[cur.items,2],
                          k27ac[cur.items,3], k27ac[cur.items,4],
                          k27ac[cur.items,5]))), col = cbbPalette[2],
          ylim = c(0,65), lwd = 8, type="o", pch=16)

    error.bars(stuff.k27me3, eyes = FALSE, sd = FALSE, bars = FALSE,
               arrow.col = cbbPalette[3], ylim = c(0,65),
               col = cbbPalette[3], add = TRUE)
    lines(colMeans((cbind(k27me3[cur.items,1], k27me3[cur.items,2],
                          k27me3[cur.items,3], k27me3[cur.items,4],
                          k27me3[cur.items,5]))), col = cbbPalette[3],
          ylim = c(0,65), lwd = 8, type="o", pch=16)

}


dev.off()
warnings()

q()
