##########
#name:          plot_feature_region_widths.r
#description:   plot widths of feature regions
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 10, 2019
##########


library(reshape2)
library(ggplot2)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

setwd(getwd())

regions.all <- read.delim(args[1], header=F, sep="\t")
regions.prom <- read.delim(args[2], header=F, sep="\t")
regions.enh <- read.delim(args[3], header=F, sep="\t")
rlength <- as.numeric(args[4])

colnames(regions.all) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
regions.all.widths = regions.all$chromEnd - regions.all$chromStart

colnames(regions.prom) <- c("chrom", "chromStart", "chromEnd", "name", 
                            "summit")
regions.prom.widths = regions.prom$chromEnd - regions.prom$chromStart

colnames(regions.enh) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
regions.enh.widths = regions.enh$chromEnd - regions.enh$chromStart

allwidths = list("promoters"= regions.prom.widths, 
                 "enhancers"=regions.enh.widths, "all"=regions.all.widths) 

allwidths.melted = melt(allwidths)
allwidths.melted$L1 = as.factor(allwidths.melted$L1)
allwidths.melted$L1 = factor(allwidths.melted$L1, 
                             levels=c(names(allwidths)), ordered=TRUE)

medians=aggregate(value ~ L1, allwidths.melted, median)
print(medians)

lengths=aggregate(value ~ L1, allwidths.melted, length)
print(lengths)

pdf(paste0("feature_regions_", rlength, "_widths.pdf"), width=8, height=6)

ggplot(allwidths.melted, aes(x=L1, y=value, fill=L1)) + 
    geom_violin(trim=TRUE, scale="count") + 
    labs(x="feature regions",y ="width in bp") + 
    theme(panel.grid=element_blank(), panel.background=element_blank(), 
          panel.border=element_rect(fill=NA), legend.position="none") + 
    scale_fill_manual(values=c("red", "darkgreen", "orange")) +
    geom_boxplot(width=0.1, outlier.shape=NA) +
    geom_text(data=medians, aes(label=value), hjust=-0.55)+
    geom_text(data=lengths, aes(y=0,label=paste0("n=",value)), 
              vjust=0, angle=45) +
    ggtitle(paste0("Feature regions widths")) +
    theme(plot.title=element_text(hjust=0.5))

dev.off()
warnings()


q()

