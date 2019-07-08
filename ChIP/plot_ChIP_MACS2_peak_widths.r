library(reshape2)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
out.dir <- args[1]
hist <- args[2]

print(hist)

setwd(getwd())

peaksD0 <- read.delim(paste0(out.dir, "/ChIP/", hist, "/D0/MACS2_peaks/bam_pooled/D0_pooled_q005_peaks.narrowPeak"), header=F, sep="\t")

peaksD2 <- read.delim(paste0(out.dir, "/ChIP/", hist, "/D2/MACS2_peaks/bam_pooled/D2_pooled_q005_peaks.narrowPeak"), header=F, sep="\t")

peaksD5 <- read.delim(paste0(out.dir, "/ChIP/", hist, "/D5/MACS2_peaks/bam_pooled/D5_pooled_q005_peaks.narrowPeak"), header=F, sep="\t")

peaksD7 <- read.delim(paste0(out.dir, "/ChIP/", hist, "/D7/MACS2_peaks/bam_pooled/D7_pooled_q005_peaks.narrowPeak"), header=F, sep="\t")

peaksD10 <- read.delim(paste0(out.dir, "/ChIP/", hist, "/D10/MACS2_peaks/bam_pooled/D10_pooled_q005_peaks.narrowPeak"), header=F, sep="\t")


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

allwidths = list("D0"= peaksD0.widths, "D2"=peaksD2.widths, "D5"=peaksD5.widths, "D7"= peaksD7.widths, "D10"=peaksD10.widths)

allwidths.melted = melt(allwidths)

allwidths.melted$L1 = as.factor(allwidths.melted$L1)
allwidths.melted$L1 = factor(allwidths.melted$L1, levels=c(names(allwidths)), ordered=TRUE)

medians=aggregate(value ~ L1, allwidths.melted, median)
print(medians)

lengths=aggregate(value ~ L1, allwidths.melted, length)
print(lengths)


pdf(paste0(hist,"_MACS2_peaks.pdf"), width=8, height=6)

ggplot(allwidths.melted, aes(x=L1, y=value, fill=L1)) + geom_violin(trim=TRUE, scale="count") + stat_summary(fun.y=median, geom="point", size=2, color="black") + labs(x=paste0(hist," peaks"), y ="lengths in bp") + theme(panel.grid=element_blank(), panel.background=element_blank(), panel.border=element_rect(fill=NA), legend.position="none", axis.text.x = element_text(angle = 0, size=15, face="bold"), axis.text.y=element_text(size=15, face="bold"), axis.ticks.y=element_line(size=2), axis.ticks.x=element_line(size=2), axis.title.x=element_text(siz=18), axis.title.y=element_text(size=18)) + scale_fill_manual(values=c("gray", "gray", "gray","gray", "gray")) +
##ggtitle("Widths distribution of ATAC-seq peaks and merged regions") +
geom_boxplot(width=0.1, outlier.shape=NA) +
geom_text(data=medians, aes(label=value), hjust=-0.55, size=4)+
geom_text(data=lengths, aes(y=0,label=paste0("n=",value)), size=4, vjust=-10,hjust=-2,angle=45)

dev.off()
warnings()


q()


