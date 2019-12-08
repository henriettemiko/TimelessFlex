##########
#name:  	plot_ATAC_MACS2_peak_widths_D8.r
#description:   plot widths of MACS2 peaks
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 27, 2019
##########


library(reshape2)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
out.dir <- args[1]
quality.dir <- args[2]

setwd(getwd())

peaksD0 <- read.delim(paste0(out.dir, 
                             "/ATAC/D0/MACS2_peaks/D0_pooled_p005_peaks_",
                             "IDR_summitpos_final.narrowPeak"), 
                      header=F, sep="\t")

peaksD2 <- read.delim(paste0(out.dir, 
                             "/ATAC/D2/MACS2_peaks/D2_pooled_p005_peaks_",
                             "IDR_summitpos_final.narrowPeak"), 
                      header=F, sep="\t")

peaksD5 <- read.delim(paste0(out.dir, 
                             "/ATAC/D5/MACS2_peaks/D5_pooled_p005_peaks_",
                             "IDR_summitpos_final.narrowPeak"), 
                      header=F, sep="\t")

peaksD7 <- read.delim(paste0(out.dir, 
                             "/ATAC/D7/MACS2_peaks/D7_pooled_p005_peaks_",
                             "IDR_summitpos_final.narrowPeak"), 
                      header=F, sep="\t")

peaksD10 <- read.delim(paste0(out.dir, 
                              "/ATAC/D10/MACS2_peaks/D10_pooled_p005_",
                              "peaks_IDR_summitpos_final.narrowPeak"), 
                       header=F, sep="\t")

peaksD8 <- read.delim(paste0(out.dir, 
                              "/ATAC/D8/MACS2_peaks/D8_pooled_p005_",
                              "peaks_IDR_summitpos_final.narrowPeak"), 
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


colnames(peaksD8) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaksD8.widths = peaksD8$chromEnd - peaksD8$chromStart

allwidths = list("D0"= peaksD0.widths, "D2"=peaksD2.widths, 
                 "D5"=peaksD5.widths, "D7"= peaksD7.widths, 
                 "D10"=peaksD10.widths, "D8"=peaksD8.widths)

allwidths.melted = melt(allwidths)

allwidths.melted$L1 = as.factor(allwidths.melted$L1)
allwidths.melted$L1 = factor(allwidths.melted$L1, 
                             levels=c(names(allwidths)), ordered=TRUE)

medians=aggregate(value ~ L1, allwidths.melted, median)

lengths=aggregate(value ~ L1, allwidths.melted, length)


pdf(paste0(quality.dir, "/ATAC_MACS2_peaks_widths_D8.pdf"), height=6, width=8)

ggplot(allwidths.melted, aes(x=L1, y=value, fill=L1)) + 
geom_violin(trim=TRUE, scale="count") + 
labs(x="time point", y ="width in bp") + 
theme(panel.grid=element_blank(), panel.background=element_blank(), 
      panel.border=element_rect(fill=NA), legend.position="none") + 
scale_fill_manual(values=c("#e8f6f3", "#abd5ee", "#63acd8","#55b7bf", 
                           "#12a8b4","#256d73")) + 
geom_boxplot(width=0.1, outlier.shape=NA) + 
geom_text(data=medians, aes(label=value), hjust=-0.55) + 
geom_text(data=lengths, aes(y=0,label=paste0("n=",value)), vjust=-20, 
          hjust=-4,angle=45) + 
ggtitle("ATAC-seq peak widths") + 
theme(plot.title=element_text(hjust=0.5))

dev.off()
warnings()


q()

