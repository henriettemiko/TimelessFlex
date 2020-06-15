##########
#name:		plot_ChIP_JAMM_peak_widths.r
#description:   plot widths of JAMM peaks
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          March 28, 2020
##########


library(reshape2)
library(ggplot2)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
out.dir <- args[1]
hist <- args[2]
quality.dir <- args[3]

print(hist)

setwd(getwd())

peaksCMP <- read.delim(paste0(out.dir, "/ChIP/", hist, "/CMP/JAMM_peaks/" , 
                             "peaks/filtered.peaks.narrowPeak"), 
                      header=F, sep="\t")

peaksMEP <- read.delim(paste0(out.dir, "/ChIP/", hist, "/MEP/JAMM_peaks/", 
                             "peaks/filtered.peaks.narrowPeak"), 
                      header=F, sep="\t")

peaksEryA <- read.delim(paste0(out.dir, "/ChIP/", hist, "/EryA/JAMM_peaks/", 
                             "peaks/filtered.peaks.narrowPeak"), 
                      header=F, sep="\t")

peaksGMP <- read.delim(paste0(out.dir, "/ChIP/", hist, "/GMP/JAMM_peaks/", 
                             "peaks/filtered.peaks.narrowPeak"), 
                      header=F, sep="\t")

peaksGranu <- read.delim(paste0(out.dir, "/ChIP/", hist, "/Granu/JAMM_peaks/", 
                             "peaks/filtered.peaks.narrowPeak"), 
                      header=F, sep="\t")

peaksMono <- read.delim(paste0(out.dir, "/ChIP/", hist, "/Mono/JAMM_peaks/", 
                             "peaks/filtered.peaks.narrowPeak"), 
                      header=F, sep="\t")


colnames(peaksCMP) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaksCMP.widths = peaksCMP$chromEnd - peaksCMP$chromStart

colnames(peaksMEP) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaksMEP.widths = peaksMEP$chromEnd - peaksMEP$chromStart

colnames(peaksEryA) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaksEryA.widths = peaksEryA$chromEnd - peaksEryA$chromStart

colnames(peaksGMP) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaksGMP.widths = peaksGMP$chromEnd - peaksGMP$chromStart

colnames(peaksGranu) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaksGranu.widths = peaksGranu$chromEnd - peaksGranu$chromStart

colnames(peaksMono) <- c("chrom", "chromStart", "chromEnd", "name", "summit")
peaksMono.widths = peaksMono$chromEnd - peaksMono$chromStart

allwidths = list("CMP"= peaksCMP.widths, "MEP"=peaksMEP.widths,
                 "EryA"=peaksEryA.widths, "GMP"= peaksGMP.widths,
                 "Granu"=peaksGranu.widths, "Mono"= peaksMono.widths)


allwidths.melted = melt(allwidths)

allwidths.melted$L1 = as.factor(allwidths.melted$L1)
allwidths.melted$L1 = factor(allwidths.melted$L1, 
                             levels=c(names(allwidths)), ordered=TRUE)

medians=aggregate(value ~ L1, allwidths.melted, median)
lengths=aggregate(value ~ L1, allwidths.melted, length)

pdf(paste0(quality.dir, "/", hist,"_JAMM_peaks_widths.pdf"), width=8, height=6)

ggplot(allwidths.melted, aes(x=L1, y=value, fill=L1)) + 
geom_violin(trim=TRUE, scale="count") + 
labs(x="time point",y ="width in bp") + 
theme(panel.grid=element_blank(), panel.background=element_blank(), 
      panel.border=element_rect(fill=NA), legend.position="none") + 
scale_fill_manual(values=c("#e8f6f3", "#abd5ee", "#63acd8","#55b7bf", 
                           "#12a8b4", "#3a5e61")) + 
geom_boxplot(width=0.1, outlier.shape=NA) + 
geom_text(data=medians, aes(label=value), hjust=-0.55) + 
geom_text(data=lengths, aes(y=0,label=paste0("n=",value)), vjust=-20, 
          hjust=-4,angle=45) + 
ggtitle(paste0(hist, " peak widths")) + 
theme(plot.title=element_text(hjust=0.5))


dev.off()
warnings()

q()


