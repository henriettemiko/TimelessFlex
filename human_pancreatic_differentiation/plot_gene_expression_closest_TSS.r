##########
#name:          plot_gene_expression_closest_TSS.r
#description:   plots gene expression of gene with closest TSS
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          January 25, 2020
##########


library(ggplot2)
library(reshape2)


fpkms = read.table("join_regions3_closestTSS_gene.txt", header=F, sep = "\t")


head(fpkms)
numTimePoints=4
numReplicates=3
numClusters=1
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

numClusters=1


pdf("FPKM_geomav_closest_TSS_region3.pdf", height=6, width=4) 
#for (i in 1:numClusters) {
  
  #cur.items = which(geomav.fpkm.reordered[[1]] == i)
  cur.scores <- geomav.fpkm.reordered[,2:5]
  
  boxplot(log(cur.scores+1), main=paste0("closest TSS to cluster 3"), 
          names = c("D0", "D2", "D5", "D10"), outline = F, 
          ylab="log(geomav gene FPKM+1)", ylim=c(0,8), col="lightskyblue1")
  
#}
  

dev.off() 
   




fpkms = read.table("join_regions7_closestTSS_gene.txt", header=F, sep = "\t")


head(fpkms)
numTimePoints=4
numReplicates=3
numClusters=1
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

numClusters=1


pdf("FPKM_geomav_closest_TSS_region7.pdf", height=6, width=4) 
#for (i in 1:numClusters) {

#cur.items = which(geomav.fpkm.reordered[[1]] == i)
cur.scores <- geomav.fpkm.reordered[,2:5]

boxplot(log(cur.scores+1), main=paste0("closest TSS to cluster 7"), 
        names = c("D0", "D2", "D5", "D10"), outline = F, 
        ylab="log(geomav gene FPKM+1)", ylim=c(0,8), col="lightskyblue1")

#}


dev.off() 

