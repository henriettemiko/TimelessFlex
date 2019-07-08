##########
#name:		plot_ChIP_JAMM_peak_widths.r
#description:   plot widths of JAMM peaks
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


library(reshape2)
library(ggplot2)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
out.dir <- args[1]
hist <- args[2]

print(hist)

setwd(getwd())

peaksD0 <- read.delim(paste0(out.dir, "/ChIP/", hist, "/D0/JAMM_peaks/" , \
                             "peaks/filtered.peaks.narrowPeak"), \
                      header=F, sep="\t")

peaksD2 <- read.delim(paste0(out.dir, "/ChIP/", hist, "/D2/JAMM_peaks/", \
                             "peaks/filtered.peaks.narrowPeak"), \
                      header=F, sep="\t")

peaksD5 <- read.delim(paste0(out.dir, "/ChIP/", hist, "/D5/JAMM_peaks/", \
                             "peaks/filtered.peaks.narrowPeak"), \
                      header=F, sep="\t")

peaksD7 <- read.delim(paste0(out.dir, "/ChIP/", hist, "/D7/JAMM_peaks/", \
                             "peaks/filtered.peaks.narrowPeak"), 
                      header=F, sep="\t")

peaksD10 <- read.delim(paste0(out.dir, "/ChIP/", hist, "/D10/JAMM_peaks/", \
                              "peaks/filtered.peaks.narrowPeak"), 
                       header=F, sep="\t")


colnames(peaksD0) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaksD0.widths = peaksD0$chromEnd - peaksD0$chromStart

colnames(peaksD2) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaksD2.widths = peaksD2$chromEnd - peaksD2$chromStart

colnames(peaksD5) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaksD5.widths = peaksD5$chromEnd - peaksD5$chromStart

colnames(peaksD7) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaksD7.widths = peaksD7$chromEnd - peaksD7$chromStart

colnames(peaksD10) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaksD10.widths = peaksD10$chromEnd - peaksD10$chromStart

allwidths = list("D0"= peaksD0.widths, "D2"=peaksD2.widths, \
                 "D5"=peaksD5.widths, "D7"= peaksD7.widths, \
                 "D10"=peaksD10.widths)

allwidths.melted = melt(allwidths)

allwidths.melted$L1 = as.factor(allwidths.melted$L1)
allwidths.melted$L1 = factor(allwidths.melted$L1, \
                             levels=c(names(allwidths)), ordered=TRUE)

medians=aggregate(value ~ L1, allwidths.melted, median)
lengths=aggregate(value ~ L1, allwidths.melted, length)

pdf(paste0(hist,"_JAMM_peaks_widths.pdf"), width=8, height=6)

ggplot(allwidths.melted, aes(x=L1, y=value, fill=L1)) + \
geom_violin(trim=TRUE, scale="count") + \
labs(x="time point",y ="width in bp") + \
theme(panel.grid=element_blank(), panel.background=element_blank(), \
      panel.border=element_rect(fill=NA), legend.position="none") + \
scale_fill_manual(values=c("#e8f6f3", "#abd5ee", "#63acd8","#55b7bf", \
                           "#12a8b4")) + \
geom_boxplot(width=0.1, outlier.shape=NA) + \
geom_text(data=medians, aes(label=value), hjust=-0.55) + \
geom_text(data=lengths, aes(y=0,label=paste0("n=",value)), vjust=-20, \
          hjust=-4,angle=45) + \
ggtitle(paste0(hist, " peak widths")) + \
theme(plot.title=element_text(hjust=0.5))


dev.off()
warnings()

q()


