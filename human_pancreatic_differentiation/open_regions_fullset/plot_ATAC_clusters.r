
##########
#name:          plot_ATAC_clusters.r
#description:   plots ATAC signal for clusters for fullset
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 11, 2019
##########


library(ggplot2)
library(reshape2)
library(plyr)
library(limma)

args = commandArgs(trailingOnly=TRUE)

curdir <- getwd()
setwd(curdir)

overlaps.prom = read.table(args[1], header=F, sep = "\t")
numTimePoints=as.numeric(args[2])
numReplicates=as.numeric(args[3])
numClusters=as.numeric(args[4])
model.dir=getwd()
numcluster.dir=getwd()

print(numReplicates)
print(numTimePoints)
print(numClusters)
print(head(overlaps.prom))


pdf(paste0(numcluster.dir,"/ATAC_cutsites_", numClusters, "_normalized.pdf"), 
    height=6, width=4)

#quantile normalization over time (rows are regions, cols are time points)
#this replaces normalization for number of cut sites
overlaps.prom.normalized = normalizeQuantiles(overlaps.prom[,2:6])

for (i in 1:numClusters) {

    cur.items.prom = which(overlaps.prom[[1]] == i)
    print(length(cur.items.prom))

    cur.overlaps.prom <- overlaps.prom.normalized[cur.items.prom,]
    print(dim(cur.overlaps.prom))

    #divide by 2 because 2 replicates were summed up
    boxplot(1000*(cur.overlaps.prom/2), 
            main=paste0("Cluster ", i, " (n=" , length(cur.items.prom), ")"), 
            ylab="normalized mean ATAC-seq cut sites",xlab="time", 
            ylim=c(0,25), col=rep("orange",5), 
            names = c("D0","D2","D5","D7","D10"),  outline = F)

}

dev.off()
warnings()


#########
#plot 2


pdf(paste0(numcluster.dir,"/ATAC_merged_peaks_", numClusters, 
           "_normalized_sep.pdf"), height=6, width=4)

for (i in 1:numClusters) {

    prom.peaks.numbers = read.table(paste0(numcluster.dir,
                                           "/merged_peaks_numbers_regions_",
                                           numClusters, "_cluster",i,".txt"), 
                                    header=F)
    print(prom.peaks.numbers)

    normvalues <- read.table(paste0("normalization_values_", numClusters, 
                                    ".txt"))

    normD0 <- normvalues[1,1]
    normD10 <- normvalues[2,1]
    normD2 <- normvalues[3,1]
    normD5 <- normvalues[4,1]
    normD7 <- normvalues[5,1]

    prom.peaks.numbers.reordered=c(prom.peaks.numbers[
                                   prom.peaks.numbers$V2=="D0",]$V1/normD0,
    prom.peaks.numbers[prom.peaks.numbers$V2=="D2",]$V1/normD2,
    prom.peaks.numbers[prom.peaks.numbers$V2=="D5",]$V1/normD5,
    prom.peaks.numbers[prom.peaks.numbers$V2=="D7",]$V1/normD7,
    prom.peaks.numbers[prom.peaks.numbers$V2=="D10",]$V1/normD10)
    print(prom.peaks.numbers.reordered)


    #for 12 cluster 0.25
    #for 8 cluster 0.4
    barplot(prom.peaks.numbers.reordered, 
            names=c("D0","D2","D5", "D7", "D10"), 
            main=paste0("Cluster ",i), xlab="time", 
            ylab="percentage of ATAC-seq peaks", ylim=c(0,0.3), col="orange")

    #for all regions together how many peaks came from which time point for 
    #merging (number of peaks/all peaks per time point)

}
dev.off()
warnings()


############
#plot 5

pdf(paste0(numcluster.dir,"/ATAC_merged_peaks_", numClusters, 
           "_fromwhichtimepoint.pdf"), height=6, width=4)

for (i in 1:numClusters) {

    #how many peaks from which time point were used for merged regions
    merged.peaks.prom = read.table(paste0(numcluster.dir,
                                          "/merged_peaks_regions_",
                                          numClusters,"_cluster",i,
                                          "_all.txt"), header=F)
    colnames(merged.peaks.prom)=c("D0", "D2", "D5", "D7", "D10")

    print(head(merged.peaks.prom))
    print(min(as.matrix(merged.peaks.prom)))
    print(max(as.matrix(merged.peaks.prom)))

    #for each region how many peaks from which time point were used for merging
    boxplot(merged.peaks.prom, names = c("D0", "D2", "D5", "D7", "D10"),  
            outline = T, col="orange", ylim=c(0,3), 
            ylab="number of merged peaks per region", yaxt='n', 
            main=paste0("Cluster ",i), xlab="time")
    axis(side=2, labels=c(0,1,2,3,4), at=c(0,1,2,3,4))


}

dev.off()
warnings()


q()

