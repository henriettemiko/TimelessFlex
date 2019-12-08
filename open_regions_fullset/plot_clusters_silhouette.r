##########
#name:          plot_clusters_silhouette.sh
#description:   plot silhouette index for fullset
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


l = read.table(paste0(numcluster.dir, "/classes-", numClusters, ".txt"))

#silhoutte index on fold changes
data1 = read.table(paste0(signalgenerator.dir, "/allFold.txt"))

data1.counts = data1[,-(1:11)]

sil1 <- silhouette(l[[1]], dist(data1.counts))
fviz_silhouette(sil1)

ggsave(paste0(type, "_", numClusters, "_silhouette_fold.pdf"), height=6, width=8)


#silhoutte index on counts
data2 = read.table(paste0(signalgenerator.dir, "/allCountsNorm.txt"))

data2.counts = data2[,-(1:11)]

sil2 <- silhouette(l[[1]], dist(data2.counts))
fviz_silhouette(sil2)

ggsave(paste0(type, "_", numClusters, "_silhouette_counts.pdf"), height=6, width=8)

warnings()
q()

