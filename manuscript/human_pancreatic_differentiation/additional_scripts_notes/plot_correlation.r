##########
#name:          plot_correlation.r
#description:   compute Spearman correlation for prom-enh initialization clusters
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          June 15, 2020
##########


library(ggplot2)
library(reshape2)
library(plyr)
library(limma)
library(psych)


#files needed:
#/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/\
#timeless_init_prom-enh/10/allCountsNorm_10classes.txt #K27ac
#all_join_leftgene_regions_fpkms.txt #FPKMs for feature regions
#/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/\
#Signal_Generator_init_prom-enh/allCountsNorm.txt
#/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/\
#timeless_init_prom-enh/2_30/classes-10.txt
#/fast/AG_Ohler/henriette/PANCREAS_final/output_hg19/open_regions_subset/open_regions_pairs/\
#timeless_init_prom-enh/10/enh_regions_10_cutsites_all.bed #ATAC


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

#for cluster 3 and 7 and noise cluster 10

for (i in c(3, 7, 10)){

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



#read in correlation values for all 3 clusters and plot

#K27ac and ATAC correlation

corrvals_c3=data.frame(read.table("Spearman_correlation_K27ac_ATAC_cluster3.txt", header=F))
corrvals_c7=data.frame(read.table("Spearman_correlation_K27ac_ATAC_cluster7.txt", header=F))
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



#K27ac and RNA correlation

corrvalsprom_c3=data.frame(read.table("Spearman_correlation_K27ac_RNA_cluster3.txt", header=F))
corrvalsprom_c7=data.frame(read.table("Spearman_correlation_K27ac_RNA_cluster7.txt", header=F))
corrvalsprom_c10=data.frame(read.table("Spearman_correlation_K27ac_RNA_cluster10.txt", header=F))


pdf("Correlation_enh_K27ac_RNA_Spearman_boxplot.pdf", width=8, height=6)
par(mfrow=c(1,3))

par(cex.main=2)
par(cex.lab=1.5)
par(cex.axis=1.5)

boxplot(corrvalsprom_c7$V1, outline=T, horizontal=F,axes=T, staplewex=1, 
        main="Cluster 7", ylab = "Spearman correlation coefficient", col="gray")
text(y = boxplot.stats(corrvalsprom_c7$V1)$stats[c(3)], 
     labels = boxplot.stats(corrvalsprom_c7$V1)$stats[c(3)], x = 1.26)

boxplot(corrvalsprom_c3$V1, outline=T, horizontal=F,axes=T, staplewex=1, 
        main="Cluster 3", ylab = "Spearman correlation coefficient", col="gray")
text(y = boxplot.stats(corrvalsprom_c3$V1)$stats[c(3)], 
     labels = boxplot.stats(corrvalsprom_c3$V1)$stats[c(3)], x = 1.26)

boxplot(corrvalsprom_c10$V1, outline=T, horizontal=F,axes=T, staplewex=1, 
        main="Cluster 10", ylab = "Spearman correlation coefficient", col="gray")
text(y = boxplot.stats(corrvalsprom_c10$V1)$stats[c(3)], 
     labels = boxplot.stats(corrvalsprom_c10$V1)$stats[c(3)], x = 1.26)


dev.off()

