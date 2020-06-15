##########
#name:          plot_closest_TSS_dist.r
#description:   plot distance to closest TSS
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 10, 2019
##########


library(reshape2)
library(ggplot2)
library(RColorBrewer)

setwd(getwd())

distances.plus <- read.delim("peaks_closestTSS_plus_dist.txt", 
                             header=F, sep="\t")
dist.plus <- distances.plus$V1
distances.minus <- read.delim("peaks_closestTSS_minus_dist.txt", 
                              header=F, sep="\t")
dist.minus <- distances.minus$V1


pdf("open_regions_closestTSS_dist.pdf", width=8, height=6)

par(mfrow=c(1,2),oma=c(0,0,2,0))

boxplot(dist.minus, outline=F, horizontal=F, axes=T, staplewex=1, 
        main="downstream edge", ylab = "distance in bp", col="gray")
text(y = boxplot.stats(dist.minus)$stats[c(1,3,5)], 
     labels = boxplot.stats(dist.minus)$stats[c(1,3,5)], x = 1.3)
text(y = boxplot.stats(dist.minus)$stats[c(2,4)], 
     labels = boxplot.stats(dist.minus)$stats[c(2,4)], x = 0.7)

boxplot(dist.plus, outline=F, horizontal=F, axes=T, staplewex=1, 
        main="upstream edge", ylab = "distance in bp", col="gray")
text(y = boxplot.stats(dist.plus)$stats[c(1,3,5)], 
     labels = boxplot.stats(dist.plus)$stats[c(1,3,5)], x = 1.3)
text(y = boxplot.stats(dist.plus)$stats[c(2,4)], 
     labels = boxplot.stats(dist.plus)$stats[c(2,4)], x = 0.7)
mtext("Distance from edges of open regions to closest TSSs", 
      outer=TRUE, cex=1.5)

dev.off()


q()

