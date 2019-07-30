
##########
#name:          plot_ATAC_clusters_pairs.r
#description:   plots ATAC signal for clusters for pairs
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 24, 2019
##########


library(ggplot2)
library(reshape2)
library(plyr)
library(limma)

args = commandArgs(trailingOnly=TRUE)

curdir <- getwd()
numcluster.dir=curdir
setwd(curdir)

overlaps.prom = read.table(args[1], header=F, sep = "\t")
overlaps.enh = read.table(args[2], header=F, sep = "\t")
numTimePoints=as.numeric(args[3])
numReplicates=as.numeric(args[4])
numClusters=as.numeric(args[5])
type=args[6]


print(head(overlaps.prom))
print(head(overlaps.enh))

    if (type=="init_prom-enh" || type=="multi_prom-enh"){
        left="PROMOTER"
        right="ENHANCER"
    } else if (type=="init_prom-prom" || type=="multi_prom-prom"){
        left="PROMOTER"
        right="PROMOTER"
    } else if (type=="init_enh-enh" || type=="multi_enh-enh"){
        left="ENHANCER"
        right="ENHANCER"
    }



pdf(paste0(numcluster.dir,"/ATAC_cutsites_", numClusters, "_normalized.pdf"), 
    height=6, width=8)

par(oma=c(0,0,2,0), xpd=TRUE)

par(mfrow=c(1,2))

#quantile normalization over time (rows are regions, cols are time points)
#this replaces normalization for number of cut sites
overlaps.prom.normalized = normalizeQuantiles(overlaps.prom[,2:5])
overlaps.enh.normalized = normalizeQuantiles(overlaps.enh[,2:5])

for (i in 1:numClusters) {

    cur.items.prom = which(overlaps.prom[[1]] == i)
    print(length(cur.items.prom))

    cur.overlaps.prom <- overlaps.prom.normalized[cur.items.prom,]
    print(dim(cur.overlaps.prom))


    cur.items.enh = which(overlaps.enh[[1]] == i)
    print(length(cur.items.enh))

    cur.overlaps.enh <- overlaps.enh.normalized[cur.items.enh,]
    print(dim(cur.overlaps.enh))


    #divide by 2 because 2 replicates were summed up
    boxplot(1000*(cur.overlaps.prom/2), 
            main=paste0(left), 
            ylab="normalized mean ATAC-seq cut sites",xlab="time", 
            ylim=c(0,25), col=rep("orange",5), 
            names = c("D0","D2","D5","D10"),  outline = F)

    boxplot(1000*(cur.overlaps.enh/2), 
            main=paste0(right), 
            ylab="normalized mean ATAC-seq cut sites",xlab="time", 
            ylim=c(0,25), col=rep("orange",5), 
            names = c("D0","D2","D5","D10"),  outline = F)

    mtext(paste0("Cluster ", i), 
          outer=TRUE, cex=1.5, font=2)


}

dev.off()
warnings()


#########
#plot 2


pdf(paste0(numcluster.dir,"/ATAC_merged_peaks_", numClusters, 
           "_normalized_sep.pdf"), height=6, width=8)

par(oma=c(0,0,2,0), xpd=TRUE)

par(mfrow=c(1,2))

for (i in 1:numClusters) {

    prom.peaks.numbers=read.table(paste0(numcluster.dir,
                                         "/prom_merged_peaks_numbers_regions_",
                                         numClusters, "_cluster",i,".txt"), 
                                  header=F)
    print(prom.peaks.numbers)

    enh.peaks.numbers=read.table(paste0(numcluster.dir,
                                        "/enh_merged_peaks_numbers_regions_",
                                        numClusters, "_cluster",i,".txt"), 
                                 header=F)
    print(enh.peaks.numbers)

    prom.normvalues <- read.table(paste0("prom_normalization_values_", 
                                         numClusters, ".txt"))

    prom.normD0 <- prom.normvalues[1,1]
    prom.normD10 <- prom.normvalues[2,1]
    prom.normD2 <- prom.normvalues[3,1]
    prom.normD5 <- prom.normvalues[4,1]

    prom.peaks.numbers.reordered=c(prom.peaks.numbers[
                                   prom.peaks.numbers$V2=="D0",]$V1/prom.normD0,
    prom.peaks.numbers[prom.peaks.numbers$V2=="D2",]$V1/prom.normD2,
    prom.peaks.numbers[prom.peaks.numbers$V2=="D5",]$V1/prom.normD5,
    prom.peaks.numbers[prom.peaks.numbers$V2=="D10",]$V1/prom.normD10)
    print(prom.peaks.numbers.reordered)


    enh.normvalues <- read.table(paste0("enh_normalization_values_", 
                                        numClusters, ".txt"))

    enh.normD0 <- enh.normvalues[1,1]
    enh.normD10 <- enh.normvalues[2,1]
    enh.normD2 <- enh.normvalues[3,1]
    enh.normD5 <- enh.normvalues[4,1]

    enh.peaks.numbers.reordered=c(enh.peaks.numbers[
                                   enh.peaks.numbers$V2=="D0",]$V1/enh.normD0,
    enh.peaks.numbers[enh.peaks.numbers$V2=="D2",]$V1/enh.normD2,
    enh.peaks.numbers[enh.peaks.numbers$V2=="D5",]$V1/enh.normD5,
    enh.peaks.numbers[enh.peaks.numbers$V2=="D10",]$V1/enh.normD10)
    print(enh.peaks.numbers.reordered)



    #for 12 cluster 0.25
    #for 8 cluster 0.4
    barplot(prom.peaks.numbers.reordered, 
            names=c("D0", "D2", "D5", "D10"), 
            main=paste0(left), xlab="time", 
            ylab="percentage of ATAC-seq peaks", ylim=c(0,0.3), col="orange")

    barplot(enh.peaks.numbers.reordered, 
            names=c("D0", "D2", "D5", "D10"), 
            main=paste0(right), xlab="time", 
            ylab="percentage of ATAC-seq peaks", ylim=c(0,0.3), col="orange")

    #for all regions together how many peaks came from which time point for 
    #merging (number of peaks/all peaks per time point)


    mtext(paste0("Cluster ", i), 
          outer=TRUE, cex=1.5, font=2)

}
dev.off()
warnings()


############
#plot 5

pdf(paste0(numcluster.dir,"/ATAC_merged_peaks_", numClusters, 
           "_fromwhichtimepoint.pdf"), height=6, width=8)

par(oma=c(0,0,2,0), xpd=TRUE)

par(mfrow=c(1,2))

for (i in 1:numClusters) {

    #how many peaks from which time point were used for merged regions
    merged.peaks.prom = read.table(paste0(numcluster.dir,
                                          "/prom_merged_peaks_regions_",
                                          numClusters,"_cluster",i,
                                          "_all.txt"), header=F)
    colnames(merged.peaks.prom)=c("D0", "D2", "D5", "D10")

    print(head(merged.peaks.prom))
    print(min(as.matrix(merged.peaks.prom)))
    print(max(as.matrix(merged.peaks.prom)))


    merged.peaks.enh = read.table(paste0(numcluster.dir,
                                          "/enh_merged_peaks_regions_",
                                          numClusters,"_cluster",i,
                                          "_all.txt"), header=F)
    colnames(merged.peaks.enh)=c("D0", "D2", "D5", "D10")

    print(head(merged.peaks.enh))
    print(min(as.matrix(merged.peaks.enh)))
    print(max(as.matrix(merged.peaks.enh)))


    #for each region how many peaks from which time point were used for merging
    boxplot(merged.peaks.prom, names = c("D0", "D2", "D5", "D10"),  
            outline = T, col="orange", ylim=c(0,3), 
            ylab="number of merged peaks per region", yaxt='n', 
            main=paste0(left), xlab="time")
    axis(side=2, labels=c(0,1,2,3,4), at=c(0,1,2,3,4))


    boxplot(merged.peaks.enh, names = c("D0", "D2", "D5", "D10"),  
            outline = T, col="orange", ylim=c(0,3), 
            ylab="number of merged peaks per region", yaxt='n', 
            main=paste0(right), xlab="time")
    axis(side=2, labels=c(0,1,2,3,4), at=c(0,1,2,3,4))

    mtext(paste0("Cluster ", i), 
          outer=TRUE, cex=1.5, font=2)

}

dev.off()
warnings()


q()

