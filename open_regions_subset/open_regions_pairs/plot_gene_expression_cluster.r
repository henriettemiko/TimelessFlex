
##########
#name:          plot_gene_expression_cluster.r
#description:   plot gene expression (FPKM from RSEM) for promoters pairs
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 12, 2019
##########


library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

curdir <- getwd()
setwd(curdir)

fpkms = read.table(args[1], header=F, sep = "\t")

numTimePoints=as.numeric(args[2])
numReplicates=as.numeric(args[3])
numClusters=as.numeric(args[4])
numcluster.dir=args[5]

numDataPoints = dim(fpkms)[1]
lastcol = dim(fpkms)[2]
mean.fpkm = data.frame(matrix(nrow=numDataPoints, ncol=numTimePoints+1))
geomav.fpkm = data.frame(matrix(nrow=numDataPoints, ncol=numTimePoints+1))

print(dim(fpkms))
#store cluster assignment in first column
mean.fpkm[,1]=fpkms[,lastcol]
geomav.fpkm[,1]=fpkms[,lastcol]

for (m in 1:numTimePoints) {
    #first column is gene id

    start=((m-1)*numReplicates)+1+1
    end=(m*numReplicates)+1
    range=seq(start, end)
    print(range)
    cur.mean <- apply(fpkms[,range], 1, function(d) (d[1]+d[2]+d[3])/3)
    cur.geomav <- apply(fpkms[,range], 1, function(d) (d[1]*d[2]*d[3])^(1/3))

    #store in table
    mean.fpkm[,m+1] = cur.mean
    geomav.fpkm[,m+1] = cur.geomav

}

#reorder to cluster,D0,D2,D5,D10 from cluster,D0,D10,D2,D5
mean.fpkm.reordered = mean.fpkm[,c(1,2,4,5,3)]
geomav.fpkm.reordered = geomav.fpkm[,c(1,2,4,5,3)]


pdf(paste0(numcluster.dir,"/FPKM_geomav_", numClusters, ".pdf"), height=6, 
    width=4) 

for (i in 1:numClusters) {

    cur.items = which(geomav.fpkm.reordered[[1]] == i)
    cur.scores <- geomav.fpkm.reordered[cur.items,2:5]

    boxplot(log(cur.scores+1), main=paste0("Cluster ", i, " PROMOTER (", 
                                           length(cur.items), " genes)"), 
            names = c("D0", "D2", "D5", "D10"), outline = F, 
            ylab="log(geomav gene FPKM+1)", ylim=c(0,8), col="lightskyblue1")

}

dev.off()
warnings()


pdf(paste0(numcluster.dir,"/FPKM_mean_", numClusters, ".pdf"), 
    height=6, width=4)

for (i in 1:numClusters) {

    cur.items = which(mean.fpkm.reordered[[1]] == i)
    cur.scores <- mean.fpkm.reordered[cur.items,2:5]

    boxplot(log(cur.scores+1), main=paste0("Cluster ", i, " PROMOTER (", 
                                           length(cur.items), " genes)"), 
            names = c("D0", "D2", "D5", "D10"), outline = F, 
            ylab="log(mean gene FPKM+1)", ylim=c(0,8), col="lightskyblue1")


}

dev.off()
warnings()


q()

