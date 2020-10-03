##########
#name:          plot_classified_open_regions_widths.r
#description:   plot widths of open regions
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          March 28, 2020
##########


library(reshape2)
library(ggplot2)
library(RColorBrewer)

setwd(getwd())

peaks.all <- read.delim("peaks_allclassified.bed", header=F, sep="\t")
colnames(peaks.all) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaks.all.widths = peaks.all$chromEnd - peaks.all$chromStart
peaks.all.dist.start = peaks.all$summit - peaks.all$chromStart
peaks.all.dist.end = peaks.all$chromEnd - peaks.all$summit

peaks.uni <- read.delim("peaks_unidirectional.bed", header=F, sep="\t")
colnames(peaks.uni) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaks.uni.widths = peaks.uni$chromEnd - peaks.uni$chromStart
peaks.uni.dist.start = peaks.uni$summit - peaks.uni$chromStart
peaks.uni.dist.end = peaks.uni$chromEnd - peaks.uni$summit

peaks.bi <- read.delim("peaks_bidirectional.bed", header=F, sep="\t")
colnames(peaks.bi) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaks.bi.widths = peaks.bi$chromEnd - peaks.bi$chromStart
peaks.bi.dist.start = peaks.bi$summit - peaks.bi$chromStart
peaks.bi.dist.end = peaks.bi$chromEnd - peaks.bi$summit

peaks.inter <- read.delim("peaks_enhancers_intergenic.bed", header=F, sep="\t")
colnames(peaks.inter) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaks.inter.widths = peaks.inter$chromEnd - peaks.inter$chromStart
peaks.inter.dist.start = peaks.inter$summit - peaks.inter$chromStart
peaks.inter.dist.end = peaks.inter$chromEnd - peaks.inter$summit

peaks.intra <- read.delim("peaks_enhancers_intragenic.bed", header=F, sep="\t")
colnames(peaks.intra) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaks.intra.widths = peaks.intra$chromEnd - peaks.intra$chromStart
peaks.intra.dist.start = peaks.intra$summit - peaks.intra$chromStart
peaks.intra.dist.end = peaks.intra$chromEnd - peaks.intra$summit


allwidths = list("unidir/div"= peaks.uni.widths, "bidir"=peaks.bi.widths, 
                 "intergenic"=peaks.inter.widths, 
                 "intragenic"=peaks.intra.widths, "all"=peaks.all.widths) 

allwidths.melted = melt(allwidths)

allwidths.melted$L1 = as.factor(allwidths.melted$L1)
allwidths.melted$L1 = factor(allwidths.melted$L1, 
                             levels=c(names(allwidths)), ordered=TRUE)

medians=aggregate(value ~ L1, allwidths.melted, median)
print(medians)

lengths=aggregate(value ~ L1, allwidths.melted, length)
print(lengths)


pdf("open_regions_widths.pdf", width=8, height=6)

ggplot(allwidths.melted, aes(x=L1, y=value, fill=L1)) + 
    geom_violin(trim=TRUE, scale="count") + 
    labs(x="open regions",y ="width in bp") + 
    theme(panel.grid=element_blank(), panel.background=element_blank(), 
          panel.border=element_rect(fill=NA), legend.position="none") + 
    scale_fill_manual(values=
                      c("red", "red", "darkgreen", "darkgreen", "orange")) +
    geom_boxplot(width=0.1, outlier.shape=NA) +
    geom_text(data=medians, aes(label=value), hjust=-0.55)+
    geom_text(data=lengths, aes(y=0,label=paste0("n=",value)), 
              vjust=0, angle=45) +
    ggtitle(paste0("Open regions widths")) +
    theme(plot.title=element_text(hjust=0.5))

dev.off()
warnings()


alldists = list("unidir"= c(peaks.uni.dist.start, peaks.uni.dist.end), 
                "bidir"=c(peaks.bi.dist.start, peaks.bi.dist.end), 
                "intergenic"=c(peaks.inter.dist.start, peaks.inter.dist.end), 
                "intragenic"=c(peaks.intra.dist.start, peaks.intra.dist.end), 
                "all"=c(peaks.all.dist.start, peaks.all.dist.end)) 

alldists.melted = melt(alldists)

alldists.melted$L1 = as.factor(alldists.melted$L1)
alldists.melted$L1 = factor(alldists.melted$L1, levels=c(names(alldists)), 
                            ordered=TRUE)

medians=aggregate(value ~ L1, alldists.melted, median)
print(medians)

lengths=aggregate(value ~ L1, alldists.melted, length)
print(lengths)


pdf("open_regions_summitdist.pdf", width=8, height=6)

ggplot(alldists.melted, aes(x=L1, y=value, fill=L1)) + 
    geom_violin(trim=TRUE, scale="count") + 
    labs(x="open regions",y ="distances from summit to edges in bp") + 
    theme(panel.grid=element_blank(), panel.background=element_blank(), 
          panel.border=element_rect(fill=NA), legend.position="none") + 
    scale_fill_manual(values=
                      c("red", "red", "darkgreen", "darkgreen", "orange")) +
    geom_boxplot(width=0.1, outlier.shape=NA) +
    geom_text(data=medians, aes(label=value), hjust=-0.55)+
    ggtitle(paste0("Distances between summit and edges of open regions")) +
    theme(plot.title=element_text(hjust=0.5))

dev.off()


pdf("open_regions_widths_boxplot.pdf", width=4, height=6)

boxplot(peaks.all.widths, outline=T, horizontal=F,axes=T, staplewex=1, 
        main="Open regions widths", ylab = "width in bp", col="gray")
text(y = boxplot.stats(peaks.all.widths)$stats[c(1,3,5)], 
     labels = boxplot.stats(peaks.all.widths)$stats[c(1,3,5)], x = 1.26)
text(y = boxplot.stats(peaks.all.widths)$stats[c(2,4)], 
     labels = boxplot.stats(peaks.all.widths)$stats[c(2,4)], x = 0.74)

dev.off()


q()

