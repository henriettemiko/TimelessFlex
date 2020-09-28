##########
#name:          plot_correlation.r
#description:   compute Spearman correlation for prom-enh initialization clusters
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          September 28, 2020
##########


library(ggplot2)
library(reshape2)
library(plyr)
library(limma)
library(psych)


#files needed:
#[hmiko@max213.mdc-berlin.net:/scratch/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_init_prom-enh/10/:
#allCountsNorm_10classes.txt #K27ac
#join_10assignments.txt #FPKM for unambigously assigned genes (not needed here)
#all_join_leftgene_regions_fpkms.txt #FPKMs for feature regions (needed here)
#/scratch/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/Signal_Generator_init_prom-enh/allCountsNorm.txt
#/scratch/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_init_prom-enh/2_30/classes-10.txt
#/scratch/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/timeless_init_prom-enh/10/enh_regions_10_cutsites_all.bed #ATAC


#K27ac signal for prom and enh (from plot_clusters_pairs.r)

score <- function(x){
  x = scale(x, scale = FALSE)
  x  = ((x-min(x))/(max(x)-min(x))) * 100
  return(x)
}


numClusters=10

l = read.table("classes-10.txt")
k = read.table("allCountsNorm.txt")

k27acprom = cbind(k[[13]], k[[14]], k[[15]], k[[16]])
k27acprom = matrix(score(as.vector(k27acprom)), ncol = 4, byrow = FALSE)

k27acenh = cbind(k[[40]], k[[41]], k[[42]], k[[43]])
k27acenh = matrix(score(as.vector(k27acenh)), ncol = 4, byrow = FALSE)

k27ac = cbind(k27acprom, k27acenh)



#RNA-seq signal for prom (from plot_expression_cluster.r for feature regions here)
numTimePoints=4
numReplicates=3

fpkms = read.table("all_join_leftgene_regions_fpkms.txt", header=F, sep = "\t")

numDataPoints = dim(fpkms)[1] #number of pairs
lastcol = dim(fpkms)[2]
geomav.fpkm = data.frame(matrix(nrow=numDataPoints, ncol=numTimePoints+1))

print(dim(fpkms))
#store cluster assignment in first column
geomav.fpkm[,1]=fpkms[,lastcol]

for (m in 1:numTimePoints) {
  #first column is gene id
  start=((m-1)*numReplicates)+1+1
  end=(m*numReplicates)+1
  range=seq(start, end)
  print(range)
  cur.geomav <- apply(fpkms[,range], 1, function(d) (d[1]*d[2]*d[3])^(1/3))
  
  #store in table
  geomav.fpkm[,m+1] = cur.geomav
}

#reorder to cluster,D0,D2,D5,D10 from cluster,D0,D10,D2,D5
geomav.fpkm.reordered = geomav.fpkm[,c(1,2,4,5,3)]



#ATAC signal for enh (from plot_ATAC_clusters_pairs.r):

overlaps.enh = read.table("enh_regions_10_cutsites_all.bed", header=F, sep = "\t")
overlaps.enh.normalized = normalizeQuantiles(overlaps.enh[,2:5])


#for all clusters, especially clusters 3 and 7 and noise cluster 10


for (i in 1:10){
  
  cur.items = which(l[[1]] == i)
  
  stuff.k27ac.prom <- cbind(k27ac[cur.items,1], k27ac[cur.items,2], 
                            k27ac[cur.items,3], k27ac[cur.items,4])
  
  k27acpromsig=stuff.k27ac.prom
  
  stuff.k27ac.enh <- cbind(k27ac[cur.items,5], k27ac[cur.items,6], 
                           k27ac[cur.items,7], k27ac[cur.items,8])
  
  k27acenhsig=stuff.k27ac.enh
  
  
  cur.items = which(geomav.fpkm.reordered[[1]] == i)
  cur.scores <- geomav.fpkm.reordered[cur.items,2:5]
  
  RNAsig=log(cur.scores+1)
  dim(RNAsig)
  
  
  cur.items.enh = which(overlaps.enh[[1]] == i)
  print(length(cur.items.enh))
  
  cur.overlaps.enh <- overlaps.enh.normalized[cur.items.enh,]
  print(dim(cur.overlaps.enh))
  
  ATACsig = 1000*(cur.overlaps.enh/2)
  
  
  
  #####Spearman correlation enh K27ac and ATAC (4-dim vectors) 
  #correlation for each feature region
  
  len=dim(k27acenhsig)[1]
  corrvals=rep(0,len)
  for (a in 1:len){
    curcorr=cor.test(k27acenhsig[a,], unlist(ATACsig[a,],use.names = F), method="spearman")
    corrvals[a]=curcorr$estimate
  }
  
  
  mean(corrvals)
  #plot(corrvals)
  #hist(corrvals)
  
  #write correlation coefficient to files
  #number of lines equals number of feature regions in cluster
  write.table(corrvals, sep = "\t",
              file = paste0("Spearman_correlation_K27ac_ATAC_cluster",i,".txt"), 
              quote = F, row.names = F, col.names = F)
  
  
  ####Spearman correlation prom K27ac and RNA (4-dim vectors)
  #correlation for each feature region
  #in case all FPKMs are zero for feature region, correlation coefficient is NA (e.g. l236 for cluster 3)
  
  len2=dim(k27acpromsig)[1]
  corrvalsprom=rep(0,len2)
  for (b in 1:len2){
    curcorrprom=cor.test(k27acpromsig[b,], unlist(RNAsig[b,],use.names = F), method="spearman")
    corrvalsprom[b]=curcorrprom$estimate
  }
  
  
  
  mean(corrvalsprom, na.rm=T)
  plot(corrvalsprom)
  hist(corrvalsprom)
  
  #write correlation coefficient to files
  #number of lines equals number of feature regions in cluster
  write.table(corrvalsprom, sep = "\t",
              file = paste0("Spearman_correlation_K27ac_RNA_cluster",i,".txt"), 
              quote = F, row.names = F, col.names = F)
  
}



#read in correlation values for all clusters and plot 3 clusters

#K27ac and ATAC correlation


corrvals_c1=data.frame(read.table("Spearman_correlation_K27ac_ATAC_cluster1.txt", header=F))
corrvals_c2=data.frame(read.table("Spearman_correlation_K27ac_ATAC_cluster2.txt", header=F))
corrvals_c3=data.frame(read.table("Spearman_correlation_K27ac_ATAC_cluster3.txt", header=F))
corrvals_c4=data.frame(read.table("Spearman_correlation_K27ac_ATAC_cluster4.txt", header=F))
corrvals_c5=data.frame(read.table("Spearman_correlation_K27ac_ATAC_cluster5.txt", header=F))
corrvals_c6=data.frame(read.table("Spearman_correlation_K27ac_ATAC_cluster6.txt", header=F))
corrvals_c7=data.frame(read.table("Spearman_correlation_K27ac_ATAC_cluster7.txt", header=F))
corrvals_c8=data.frame(read.table("Spearman_correlation_K27ac_ATAC_cluster8.txt", header=F))
corrvals_c9=data.frame(read.table("Spearman_correlation_K27ac_ATAC_cluster9.txt", header=F))
corrvals_c10=data.frame(read.table("Spearman_correlation_K27ac_ATAC_cluster10.txt", header=F))


pdf("Correlation_enh_K27ac_ATAC_Spearman_boxplot.pdf", width=8, height=6)
par(mfrow=c(1,3))

par(cex.main=2)
par(cex.lab=1.5)
par(cex.axis=1.5)

boxplot(corrvals_c7$V1, outline=T, horizontal=F,axes=T, staplewex=1, 
        main="Cluster 7", ylab = "Spearman correlation coefficient", col="gray")
text(y = boxplot.stats(corrvals_c7$V1)$stats[c(3)], 
     labels = boxplot.stats(corrvals_c7$V1)$stats[c(3)], x = 1.26)

boxplot(corrvals_c3$V1, outline=T, horizontal=F,axes=T, staplewex=1, 
        main="Cluster 3", ylab = "Spearman correlation coefficient", col="gray")
text(y = boxplot.stats(corrvals_c3$V1)$stats[c(3)], 
     labels = boxplot.stats(corrvals_c3$V1)$stats[c(3)], x = 1.26)

boxplot(corrvals_c10$V1, outline=T, horizontal=F,axes=T, staplewex=1, 
        main="Cluster 10", ylab = "Spearman correlation coefficient", col="gray")
text(y = boxplot.stats(corrvals_c10$V1)$stats[c(3)], 
     labels = boxplot.stats(corrvals_c10$V1)$stats[c(3)], x = 1.26)


dev.off()


###for each cluster get median and then plot medians


m1 = boxplot.stats(corrvals_c1$V1)$stats[c(3)]
m2 = boxplot.stats(corrvals_c2$V1)$stats[c(3)]
m3 = boxplot.stats(corrvals_c3$V1)$stats[c(3)]
m4 = boxplot.stats(corrvals_c4$V1)$stats[c(3)]
m5 = boxplot.stats(corrvals_c5$V1)$stats[c(3)]
m6 = boxplot.stats(corrvals_c6$V1)$stats[c(3)]
m7 = boxplot.stats(corrvals_c7$V1)$stats[c(3)]
m8 = boxplot.stats(corrvals_c8$V1)$stats[c(3)]
m9 = boxplot.stats(corrvals_c9$V1)$stats[c(3)]
m10 = boxplot.stats(corrvals_c10$V1)$stats[c(3)]

pdf("Correlation_enh_K27ac_ATAC_Spearman_median_all.pdf", width=8, height=6)
par(mfrow=c(1,1))

plot (1:10, c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10), xlab= "cluster", ylab = "median of Spearman correlation coefficents", ylim=c(0,1), pch=19)
axis(1, at=1:10, labels=1:10)
dev.off()

