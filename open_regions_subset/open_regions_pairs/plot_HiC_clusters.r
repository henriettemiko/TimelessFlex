

library("ggplot2")


args = commandArgs(trailingOnly=TRUE)


#if (length(args)==0) {
#  stop("One arguments must be provided: NUM_CLUSTER", call.=FALSE)
#}# else if (length(args)==1) {
##  # default output file
##  args[2] = "out.txt"
##}


interactions=read.table(args[1],header=F,stringsAsFactors=F)
interactions2=read.table(args[2],header=F,stringsAsFactors=F)
normvalues=read.table(args[3],header=F,stringsAsFactors=F)
numClusters=as.numeric(args[4])
print(numClusters)
numcluster.dir=getwd()



pdf(paste0(numcluster.dir,"/HiC_interactions_", numClusters, "_new.pdf"), height=6, width=4)


for (i in 1:numClusters) {
  
 # print(numbers[i])

  #how many Hi-C interactions do we have at each time point between the prom-enh pairs in the clusters?
  #counts how often D0, D2, D5, D10 occur in each cluster

  #interactions = read.table(paste0(numcluster.dir,"/HiC_",numClusters, "_new_cluster",i,"_newline_counts.txt"), header=F, stringsAsFactors = F)

  
  print(which(is.na(as.matrix(interactions))))
  
  print("here")
  #    615 D0
  #    723 D10
  #   1086 D2
  #    637 D5
  
 # > 0.1788618 + 0.1951220 + 0.07967480 + 0.08455285 +  0.07804878 + 0.07154472  + 0.11544715 + 0.1967480
#  [1] 1 



#new
#10
#   1236 D0
#   1419 D10
#   2405 D2
#   1253 D5



    #divide by number of total number of D0, D2, D5, D10 interactions in all clusters
  interactions.reordered=c(interactions[interactions$V2=="D0",]$V1/1236,interactions[interactions$V2=="D2",]$V1/2405,interactions[interactions$V2=="D5",]$V1/1253,interactions[interactions$V2=="D10",]$V1/1419)
  print(interactions.reordered)
 
  
# # for 8  
# #  [1] 0.07274119 0.12266773 0.07393293 0.06818182
# #  [1] 0.09111792 0.10996427 0.10594512 0.13836898
# #  [1] 0.2220521 0.1587932 0.1730183 0.1824866
# #  [1] 0.2519142 0.2346169 0.2583841 0.2613636
# #  [1] 0.09724349 0.08614530 0.10823171 0.08689840
# #  [1] 0.03062787 0.04446209 0.02057927 0.02005348
# #  [1] 0.1546708 0.1445018 0.1737805 0.1550802
# #  [1] 0.07963247 0.09884875 0.08612805 0.08756684
#   
#  # 0.07274119 + 0.09111792 + 0.2220521 + 0.2519142 + 0.09724349 + 0.03062787 + 0.1546708 + 0.07963247
#  # [1] 1
#   
#   
#   #and divide by number of pairs in cluster
#   #barplot(100*interactions.reordered/numbers[i], names=c("D0","D2","D5","D10"), main=paste0("Cluster ",i), xlab="time", ylab="Normalized percentage of Hi-C interactions", ylim=c(0,0.1))
  
  #for 8 cluster 0.3 is okay
  barplot(interactions.reordered, names=c("D0","D2","D5","D10"), main=paste0("Cluster ",i), xlab="time", ylab="percentage of Hi-C interactions at time point", ylim=c(0,0.25), col="lightskyblue")
  
 # #for 12 cluster 0.25 is okay
 # barplot(interactions.reordered, names=c("D0","D2","D5","D10"), main=paste0("Cluster ",i), xlab="time", ylab="percentage of Hi-C interactions", ylim=c(0,0.25), col="lightskyblue")
#   
#   
#  # 1306 D0
# #  1496 D10
# #  2519 D2
# #  1312 D5

}
dev.off()
warnings()


#q()

#############################################


pdf(paste0(numcluster.dir,"/HiC_interactions_", numClusters, "_forwhichtimepoint.pdf"), height=6, width=4)


for (i in 1:numClusters) {
  
  #print(numbers[i])

  #for each pair in cluster: did it had interaction at D0,D2,D5,D10? 1: yes, 0: no
  #does it have an interaction (no:0, yes:1)
  #interactions2 = read.table(paste0(numcluster.dir,"/HiC_",numClusters,"_new_cluster",i,"_all.txt"), header=F, stringsAsFactors = F)
  colnames(interactions2)=c("D0", "D2", "D5", "D10")

  print(min(as.matrix(interactions2)))
  print(max(as.matrix(interactions2)))

  #TODO: change ylim depending on max value
  
  #for each region how many peaks from which time point were used for merging
  boxplot(interactions2, names = c("D0","D2","D5","D10"),  outline = T, col="lightskyblue", ylim=c(0,1), ylab="number of Hi-C interactions per region", yaxt='n', main=paste0("Cluster ",i), xlab="time")
  axis(side=2, labels=c(0,1,2,3), at=c(0,1,2,3))
  #boxplot(interactions2, names = c("D0","D2","D5","D10"),  outline = T, main=paste0("Cluster ", i), col="lightskyblue", ylim=c(0,1), ylab="Interaction yes/no", xlab="time")
  
  print(which(is.na(as.matrix(interactions2))))
  print("here2")
  
  print(head(as.matrix(interactions2)))
  
}
dev.off()
warnings()






#there is a maximum of 1 interaction per time point







library(gplots)
library(RColorBrewer)

pdf(paste0(numcluster.dir,"/HiC_interactions_", numClusters, "_heatmap.pdf"), height=6, width=4)

for (i in 1:numClusters) {
  
  #my_palette=c("black", "orange", "red", "purple")
  my_palette=c("black", "orange")
  #col_breaks=c(0,0.55,1.01)
  #heatmap.2(as.matrix(interactions2), keysize=1.5, key.xlab = "Hi-C interactions per region", key.title  = "", key.ylab = "", main=paste0("Cluster ",i), Rowv=T,Colv=F, col=my_palette, labRow="",labCol=c("D0","D2","D5","D10"),dendrogram="none",scale="none",density.info = "none", trace="none",ylab = "regions", xlab="time",margins=c(6,2))
  
  
  
  
  heatmap.2(as.matrix(interactions2), key.xtickfun = function() {return(list(at=c(0,1), labels=c(0,1)))} ,keysize=1.5, key.xlab = "Hi-C interactions per region", key.title  = "", key.ylab = "", main=paste0("Cluster ",i), Rowv=T,Colv=F, col=my_palette, labRow="",labCol=c("D0","D2","D5","D10"),dendrogram="none",scale="none",density.info = "none", trace="none",ylab = "regions", xlab="time",margins=c(6,2))
  
  

  
  
 # my_palette=c("black", "orange", "red", "purple")
#  col_breaks=c(0,0.5,1.5,2.5,3.5)
#  heatmap.2(as.matrix(merged.peaks.enh), keysize=1.5, key.xlab = "merged peaks per region", key.title  = "", key.ylab = "",main=paste0("Cluster ",i, ": ENHANCER"), Rowv=T,Colv=F, col=my_palette, breaks=col_breaks, labRow="",labCol=c("D0","D2","D5","D10"),dendrogram="none",scale="none",density.info = "none", trace="none",ylab = "regions", xlab="time",margins=c(6,2))
  
}

dev.off()
warnings()





q()













######I THINK THIS DOES NOT MAKE SENSE

pdf(paste0(numcluster.dir,"/HiC_interactions_", numClusters, "_new3.pdf"), height=5, width=4)


for (i in 1:numClusters) {
  
  interactions2 = read.table(paste0(numcluster.dir,"/HiC_8_new_cluster",i,"_all.txt"), header=F, stringsAsFactors = F)

  colnames(interactions2)=c("D0", "D2", "D5", "D10")
  
  print(colMeans(interactions2))
#I AM NOT SURE THIS MAKES SENSE  
  #mean of does it have an interaction (no:0, yes:1)
  barplot(colMeans(interactions2), main=paste0("Cluster ", i), col="red", ylim=c(0,1), ylab="mean of interaction yes/no", xlab="time")

  
}
dev.off()
warnings()

q()



####################



#tidyverse
#  data_frame %>% filter(V2=='D0') %>% .$V1

