
##########
#name:          plot_gene_expression_cluster.r
#description:   plot gene expression for promoters from pairs prom-prom
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 24, 2019
##########


library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

curdir <- getwd()
setwd(curdir)

fpkms1 = read.table(args[1], header=F, sep = "\t")
fpkms2 = read.table(args[2], header=F, sep = "\t")

numTimePoints=as.numeric(args[3])
numReplicates=as.numeric(args[4])
numClusters=as.numeric(args[5])
numcluster.dir=args[6]

numDataPoints1 = dim(fpkms1)[1]
lastcol1 = dim(fpkms1)[2]
mean.fpkm1 = data.frame(matrix(nrow=numDataPoints1, ncol=numTimePoints+1))
geomav.fpkm1 = data.frame(matrix(nrow=numDataPoints1, ncol=numTimePoints+1))

numDataPoints2 = dim(fpkms2)[1]
lastcol2 = dim(fpkms2)[2]
mean.fpkm2 = data.frame(matrix(nrow=numDataPoints2, ncol=numTimePoints+1))
geomav.fpkm2 = data.frame(matrix(nrow=numDataPoints2, ncol=numTimePoints+1))

print(dim(fpkms1))
print(dim(fpkms2))

#store cluster assignment in first column
mean.fpkm1[,1]=fpkms1[,lastcol1]
geomav.fpkm1[,1]=fpkms1[,lastcol1]
mean.fpkm2[,1]=fpkms2[,lastcol2]
geomav.fpkm2[,1]=fpkms2[,lastcol2]


for (m in 1:numTimePoints) {
    #first column is gene id

    start=((m-1)*numReplicates)+1+1
    end=(m*numReplicates)+1
    range=seq(start, end)
    print(range)
    cur.mean1 <- apply(fpkms1[,range], 1, function(d) (d[1]+d[2]+d[3])/3)
    cur.geomav1 <- apply(fpkms1[,range], 1, function(d) (d[1]*d[2]*d[3])^(1/3))

    cur.mean2 <- apply(fpkms2[,range], 1, function(d) (d[1]+d[2]+d[3])/3)
    cur.geomav2 <- apply(fpkms2[,range], 1, function(d) (d[1]*d[2]*d[3])^(1/3))

    #store in table
    mean.fpkm1[,m+1] = cur.mean1
    geomav.fpkm1[,m+1] = cur.geomav1
    mean.fpkm2[,m+1] = cur.mean2
    geomav.fpkm2[,m+1] = cur.geomav2


}

#reorder to cluster,D0,D2,D5,D10 from cluster,D0,D10,D2,D5
mean.fpkm.reordered1 = mean.fpkm1[,c(1,2,4,5,3)]
geomav.fpkm.reordered1 = geomav.fpkm1[,c(1,2,4,5,3)]
mean.fpkm.reordered2 = mean.fpkm2[,c(1,2,4,5,3)]
geomav.fpkm.reordered2 = geomav.fpkm2[,c(1,2,4,5,3)]

pdf(paste0(numcluster.dir,"/FPKM_geomav_", numClusters, ".pdf"), height=6, 
    width=8) 

par(mfrow=c(1,2))

for (i in 1:numClusters) {

    cur.items1 = which(geomav.fpkm.reordered1[[1]] == i)
    cur.scores1 <- geomav.fpkm.reordered1[cur.items1,2:5]

    boxplot(log(cur.scores1+1), main=paste0("Cluster ", i, ": ", 
                                           length(cur.items1), " genes"), 
            names = c("D0", "D2", "D5", "D10"), outline = F, 
            ylab="log(geomav gene FPKM+1)", ylim=c(0,8), col="lightskyblue1")

    cur.items2 = which(geomav.fpkm.reordered2[[1]] == i)
    cur.scores2 <- geomav.fpkm.reordered2[cur.items2,2:5]

    boxplot(log(cur.scores2+1), main=paste0(length(cur.items2), " genes"), 
            names = c("D0", "D2", "D5", "D10"), outline = F, 
            ylab="log(geomav gene FPKM+1)", ylim=c(0,8), col="lightskyblue1")


}

dev.off()
warnings()


pdf(paste0(numcluster.dir,"/FPKM_mean_", numClusters, ".pdf"), 
    height=6, width=8)

par(mfrow=c(1,2))

for (i in 1:numClusters) {

    cur.items1 = which(mean.fpkm.reordered1[[1]] == i)
    cur.scores1 <- mean.fpkm.reordered1[cur.items1,2:5]

    boxplot(log(cur.scores1+1), main=paste0("Cluster ", i, ": ", 
                                           length(cur.items1), " genes"), 
            names = c("D0", "D2", "D5", "D10"), outline = F, 
            ylab="log(mean gene FPKM+1)", ylim=c(0,8), col="lightskyblue1")


    cur.items2 = which(mean.fpkm.reordered2[[1]] == i)
    cur.scores2 <- mean.fpkm.reordered2[cur.items2,2:5]

    boxplot(log(cur.scores2+1), main=paste0(length(cur.items2), " genes"), 
            names = c("D0", "D2", "D5", "D10"), outline = F, 
            ylab="log(mean gene FPKM+1)", ylim=c(0,8), col="lightskyblue1")

}

dev.off()
warnings()


q()

