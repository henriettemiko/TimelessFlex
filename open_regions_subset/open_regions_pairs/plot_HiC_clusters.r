
##########
#name:          plot_HiC_clusters.r
#description:   plots Hi-C signal for clusters for pairs
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 29, 2019
##########

library("ggplot2")


args = commandArgs(trailingOnly=TRUE)

numClusters=as.numeric(args[1])
print(numClusters)
numcluster.dir=getwd()


pdf(paste0(numcluster.dir,"/HiC_interactions_", numClusters, ".pdf"), 
    height=6, width=4)

for (i in 1:numClusters) {
  

    #how many Hi-C interactions do we have at each time point between the 
    #pairs in the clusters?
    #counts how often D0, D2, D5, D10 occur in each cluster
    interactions = read.table(paste0(numcluster.dir,"/HiC_",numClusters, 
                                     "_cluster",i,"_newline_counts.txt"), 
                              header=F, stringsAsFactors = F)
  
    print(which(is.na(as.matrix(interactions))))

    normvalues <- read.table(paste0("normalization_values_", numClusters, 
                                    ".txt"))
    normD0 <- normvalues[1,1]
    normD10 <- normvalues[2,1]
    normD2 <- normvalues[3,1]
    normD5 <- normvalues[4,1]

    #divide by number of total number of D0, D2, D5, D10 interactions in all 
    #clusters
    interactions.reordered=c(interactions[interactions$V2=="D0",]$V1/normD0,
                             interactions[interactions$V2=="D2",]$V1/normD2,
                             interactions[interactions$V2=="D5",]$V1/normD5,
                             interactions[interactions$V2=="D10",]$V1/normD10)


  #for 8 cluster 0.3 is okay
  barplot(interactions.reordered, names=c("D0","D2","D5","D10"), 
          main=paste0("Cluster ",i), xlab="time", 
          ylab="percentage of Hi-C interactions at time point", 
          ylim=c(0,0.25), col="lightskyblue")
  
}
dev.off()
warnings()


#############################################

pdf(paste0(numcluster.dir,"/HiC_interactions_", numClusters, 
           "_forwhichtimepoint.pdf"), height=6, width=4)


for (i in 1:numClusters) {
  
    #for each pair in cluster: did it had interaction at D0,D2,D5,D10?
    #no:0, yes:1
    interactions2 = read.table(paste0(numcluster.dir,"/HiC_",numClusters,
                                      "_cluster",i,"_all.txt"), header=F, 
                               stringsAsFactors = F)
    colnames(interactions2)=c("D0", "D2", "D5", "D10")

    print(min(as.matrix(interactions2)))
    print(max(as.matrix(interactions2)))

    #TODO: change ylim depending on max value
  
    #for each region how many peaks from which time point were used for merging
    boxplot(interactions2, names = c("D0","D2","D5","D10"),  outline = T, 
            col="lightskyblue", ylim=c(0,1), 
            ylab="number of Hi-C interactions per region", yaxt='n', 
            main=paste0("Cluster ",i), xlab="time")
    axis(side=2, labels=c(0,1,2,3), at=c(0,1,2,3))
  
    print(which(is.na(as.matrix(interactions2))))
  
    print(head(as.matrix(interactions2)))
  
}

dev.off()
warnings()


##############################################################################

library(gplots)
library(RColorBrewer)

pdf(paste0(numcluster.dir,"/HiC_interactions_", numClusters, "_heatmap.pdf"), 
    height=6, width=4)

#there is a maximum of 1 interaction per time point

for (i in 1:numClusters) {
  
    my_palette=c("black", "orange")
  
    heatmap.2(as.matrix(interactions2), 
              key.xtickfun = function(){return(list(at=c(0,1),labels=c(0,1)))},
              keysize=1.5, key.xlab = "Hi-C interactions per region", 
              key.title  = "", key.ylab = "", main=paste0("Cluster ",i), 
              Rowv=T,Colv=F, col=my_palette, labRow="",
              labCol=c("D0","D2","D5","D10"),dendrogram="none",scale="none",
              density.info = "none", trace="none",ylab = "regions", 
              xlab="time",margins=c(6,2))
  
}

dev.off()
warnings()


q()

