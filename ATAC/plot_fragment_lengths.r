##########
#name:          plot_fragment_lengths.r
#description:   plots fragment lengths
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 9, 2019
##########


library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

#curdir <- getwd()
#setwd(curdir)

counts1 = read.table(args[1], header=F)
counts2 = read.table(args[2], header=F)
time=args[3]
num.lines1=as.numeric(args[4])
num.lines2=as.numeric(args[5])
quality.dir=args[6]


#divide counts by number of reads in bedpe file
norm.counts1 <- counts1[,1]/num.lines1*10^3
norm.counts2 <- counts2[,1]/num.lines2*10^3

pdf(paste0(quality.dir,"/ATAC_",time,"_rep1_fragment_length.pdf"), 
    width=8, height=6)
plot(counts1[,2], norm.counts1, 
     main=paste0("Fragment lengths ATAC ", time, " rep1",), 
     xlab="fragment length", type="l", ylab="normalized counts")
dev.off()

pdf(paste0(quality.dir,"/ATAC_",time,"_rep2_fragment_length.pdf"), 
    width=8, height=6)
plot(counts2[,2], norm.counts2, 
     main=paste0("Fragment lengths ATAC ", time, " rep2",), 
     xlab="fragment length", type="l", ylab="normalized counts")
dev.off()


#log scale
pdf(paste0(quality.dir,"/ATAC_",time, "_rep1_fragment_length_log.pdf"), 
    width=8, height=6)
plot(counts1[,2], log(norm.counts1+1), main=paste0("Fragment lengths ATAC ", 
                                                 time, " rep1"), 
     xlab="fragment length", type="l", ylab="log(normalized counts+1)")
dev.off()

pdf(paste0(quality.dir,"/ATAC_",time, "_rep2_fragment_length_log.pdf"), 
    width=8, height=6)
plot(counts2[,2], log(norm.counts2+1), main=paste0("Fragment lengths ATAC ", 
                                                 time, " rep2"), 
     xlab="fragment length", type="l", ylab="log(normalized counts+1)")
dev.off()

q()

