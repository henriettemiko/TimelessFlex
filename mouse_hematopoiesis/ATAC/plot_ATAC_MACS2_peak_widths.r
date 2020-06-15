##########
#name:  	plot_ATAC_MACS2_peak_widths.r
#description:   plot widths of MACS2 peaks
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          April 3, 2020
##########


library(reshape2)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
out.dir <- args[1]
quality.dir <- args[2]

setwd(getwd())

peaksCMP <- read.delim(paste0(out.dir, 
                             "/ATAC/CMP/MACS2_peaks/CMP_SRR1533865_",
                             "trimmed_filtered_sorted_q001_peaks_",
                             "summitpos_final.narrowPeak"), 
                      header=F, sep="\t")

peaksMEP <- read.delim(paste0(out.dir, 
                             "/ATAC/MEP/MACS2_peaks/MEP_SRR1533868_",
                             "trimmed_filtered_sorted_q001_peaks_",
                             "summitpos_final.narrowPeak"), 
                      header=F, sep="\t")

peaksEryA <- read.delim(paste0(out.dir, 
                             "/ATAC/EryA/MACS2_peaks/EryA_SRR3001797_",
                             "trimmed_filtered_sorted_q001_peaks_",
                             "summitpos_final.narrowPeak"), 
                      header=F, sep="\t")

peaksGMP <- read.delim(paste0(out.dir, 
                             "/ATAC/GMP/MACS2_peaks/GMP_SRR1533850_",
                             "trimmed_filtered_sorted_q001_peaks_",
                             "summitpos_final.narrowPeak"), 
                      header=F, sep="\t")

peaksGranu <- read.delim(paste0(out.dir, 
                             "/ATAC/Granu/MACS2_peaks/Granu_SRR1533857_",
                             "trimmed_filtered_sorted_q001_peaks_",
                             "summitpos_final.narrowPeak"), 
                      header=F, sep="\t")

peaksMono <- read.delim(paste0(out.dir, 
                             "/ATAC/Mono/MACS2_peaks/Mono_SRR1533852_",
                             "trimmed_filtered_sorted_q001_peaks_",
                             "summitpos_final.narrowPeak"), 
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


pdf(paste0(quality.dir, "/ATAC_MACS2_peaks_widths.pdf"), height=6, width=8)

ggplot(allwidths.melted, aes(x=L1, y=value, fill=L1)) + 
geom_violin(trim=TRUE, scale="count") + 
labs(x="time point", y ="width in bp") + 
theme(panel.grid=element_blank(), panel.background=element_blank(), 
      panel.border=element_rect(fill=NA), legend.position="none") + 
scale_fill_manual(values=c("#e8f6f3", "#abd5ee", "#63acd8","#55b7bf", 
                           "#12a8b4", "#3a5e61")) + 
geom_boxplot(width=0.1, outlier.shape=NA) + 
geom_text(data=medians, aes(label=value), hjust=-0.55) + 
geom_text(data=lengths, aes(y=0,label=paste0("n=",value)), vjust=-20, 
          hjust=-4,angle=45) + 
ggtitle("ATAC-seq peak widths") + 
theme(plot.title=element_text(hjust=0.5))

dev.off()
warnings()


q()

