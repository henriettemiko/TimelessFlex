##########
#name:          plot_fragment_lengths.r
#description:   plots fragment lengths
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 9, 2019
##########


library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

curdir <- getwd()
setwd(curdir)

name=args[1]
time=args[2]
num.lines=as.numeric(args[3])
quality.dir=args[4]

print(quality.dir)
print(name)
print(num.lines)

counts = read.table(paste0(quality.dir, "/ATAC_",time, "_", name, 
                           "_fragment_lengths_counts.txt"), header=F)

#divide counts by number of reads in bedpe file
norm.counts <- counts[,1]/num.lines*10^3

pdf(paste0(quality.dir,"/ATAC_",time,"_" ,name,"_fragment_length.pdf"), 
    width=8, height=6)
plot(counts[,2], norm.counts, 
     main=paste0("Fragment lengths ATAC ", time, " ", name), 
     xlab="fragment length", type="l", ylab="normalized counts")
dev.off()


#log scale
pdf(paste0(quality.dir,"/ATAC_",time, "_",name,"_fragment_length_log.pdf"), 
    width=8, height=6)
plot(counts[,2], log(norm.counts+1), main=paste0("Fragment lengths ATAC ", 
                                                 time, " ", name), 
     xlab="fragment length", type="l", ylab="log(normalized counts+1)")
dev.off()


q()

