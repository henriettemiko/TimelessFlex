##########
#name:          plot_feature_region_widths.r
#description:   plot widths of feature regions
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 10, 2019
##########


library(reshape2)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)

setwd(getwd())

feature.dir <- args[1]
rlength <- as.numeric(args[2])
name <- args[3]
titlename <- args[4]

#read in overlaps for each region (one line is ALL overlaps with histone mark 
#for this region
#sometimes 1, sometimes more
#number of rows are number of regions that have any overlap with this histone 
#mark
#if there is more than one overlap then I take the maximum overlap for the plot

cols.K27ac <- max(count.fields(paste0(feature.dir, "/feature_regions_",
                                      rlength,"_", name, 
                                      "H3K27ac_overlaps.bed"), sep="\t"))

overlaps.K27ac <- as.matrix(read.delim(paste0(feature.dir, 
                                              "/feature_regions_",
                                              rlength,"_", name,
                                              "H3K27ac_overlaps.bed"), 
                                       fill=TRUE,header=F, sep="\t",
                                       col.names=1:cols.K27ac))
print(head(overlaps.K27ac))

cols.K27me3 <- max(count.fields(paste0(feature.dir, 
                                       "/feature_regions_",rlength,"_", 
                                       name, "H3K27me3_overlaps.bed"), 
                                sep="\t"))

overlaps.K27me3 <- as.matrix(read.delim(paste0(feature.dir, 
                                               "/feature_regions_",rlength,"_",
                                               name, "H3K27me3_overlaps.bed"), 
                                        fill=TRUE,header=F, sep="\t",
                                        col.names=1:cols.K27me3))

cols.K4me3 <- max(count.fields(paste0(feature.dir, 
                                      "/feature_regions_",rlength,"_" ,name, 
                                      "H3K4me3_overlaps.bed"), sep="\t"))

overlaps.K4me3 <- as.matrix(read.delim(paste0(feature.dir, 
                                              "/feature_regions_",rlength,"_", 
                                              name, "H3K4me3_overlaps.bed"), 
                                       fill=TRUE,header=F, sep="\t",
                                       col.names=1:cols.K4me3))

cols.K4me1 <- max(count.fields(paste0(feature.dir, 
                                      "/feature_regions_",rlength,"_", name, 
                                      "H3K4me1_overlaps.bed"), sep="\t"))

overlaps.K4me1 <- as.matrix(read.delim(paste0(feature.dir, 
                                              "/feature_regions_",rlength,"_", 
                                              name, "H3K4me1_overlaps.bed"), 
                                       fill=TRUE,header=F, sep="\t",
                                       col.names=1:cols.K4me1))

median.K27ac=apply(overlaps.K27ac, 1, max, na.rm=TRUE)
median.K27me3=apply(overlaps.K27me3, 1, max, na.rm=TRUE)
median.K4me3=apply(overlaps.K4me3, 1, max, na.rm=TRUE)
median.K4me1=apply(overlaps.K4me1, 1, max, na.rm=TRUE)

allwidths = list("H3K27ac"=median.K27ac, "H3K27me3"=median.K27me3,
                 "H3K4me1"=median.K4me1,"H3K4me3"=median.K4me3 ) 

allwidths.melted = melt(allwidths) 
allwidths.melted$L1 = as.factor(allwidths.melted$L1)
allwidths.melted$L1 = factor(allwidths.melted$L1, 
                             levels=c(names(allwidths)), ordered=TRUE)

medians=aggregate(value ~ L1, allwidths.melted, median)
print(medians)

lengths=aggregate(value ~ L1, allwidths.melted, length)
print(lengths)


histcolors <- c("#009E73", "#D50F25", "white", "gray" ) 
#me3 gray, me1 black, 27ac green 009E73, 27me3 red
#plot from dark to light colors: me1, k27me3, k27ac, me3

pdf(paste0(feature.dir,"/feature_regions_",rlength,"_", name, "overlaps.pdf"), 
    width=8, height=6)

ggplot(allwidths.melted, aes(x=L1, y=value, fill=L1)) + 
    geom_violin(trim=TRUE, scale="count") +
    labs(x=paste0(titlename, " feature regions"),y ="max overlap in bp") + 
    theme(panel.grid=element_blank(), panel.background=element_blank(), 
          panel.border=element_rect(fill=NA), legend.position="none") +
    scale_fill_manual(values=histcolors) +
    geom_boxplot(width=0.1, outlier.shape=NA) +
    geom_text(data=medians, aes(label=value), hjust=-0.55)+
    ggtitle(paste0("Maximum overlaps of ", titlename, 
                   " feature regions with histone marks")) +
    theme(plot.title=element_text(hjust=0.5))

dev.off()
warnings()


q()

