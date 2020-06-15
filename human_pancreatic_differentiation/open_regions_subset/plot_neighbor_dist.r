##########
#name:          plot_neighbor_dist.r
#description:   plot distance to neighbor
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 10, 2019
##########


library(reshape2)
library(ggplot2)
library(RColorBrewer)

setwd(getwd())

#distribution of distance between regions

neighbors <- read.delim("peaks_allclassified_neighbors.bed", 
                        header=F, sep="\t")
neighbors$dist = neighbors$V13 - neighbors$V3


pdf("open_regions_neighbors_dist.pdf", width=4, height=6)

#same chromosome! V1 is chromosome of first peak and V12 is chromosome of 
#second peak
distances <- neighbors$dist[neighbors$V1==neighbors$V12]

boxplot(distances, outline=F, horizontal=F, axes=T, staplewex=1, 
        main="Distances between open regions", ylab = "length in bp", 
        col="gray")
text(y = boxplot.stats(distances)$stats[c(1,3,5)], 
     labels = boxplot.stats(distances)$stats[c(1,3,5)], x = 1.32)
text(y = boxplot.stats(distances)$stats[c(2,4)], 
     labels = boxplot.stats(distances)$stats[c(2,4)], x = 0.7)

dev.off()


q()

