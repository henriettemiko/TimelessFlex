##########
#name:          get_normalized_FC.r
#description:   compute normalized fold changes of consecutive time points
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 10, 2019
##########

library(limma)

args = commandArgs(trailingOnly=TRUE)

curdir <- getwd()
setwd(curdir)

l = read.table(args[1])
name <- args[2]
titlename <- args[3]

print(dim(l))

#select columns for histone mark (columns = time points)
#regions have 11 columns at beginning of file
m1 = normalizeQuantiles(cbind(l[[12]], l[[13]], l[[14]],l[[15]], l[[16]], l[[17]])) 
m2 = normalizeQuantiles(cbind(l[[18]], l[[19]],l[[20]], l[[21]], l[[22]], l[[23]])) 
m3 = normalizeQuantiles(cbind(l[[24]], l[[25]], l[[26]],l[[27]], l[[28]], l[[29]])) 
m4 = normalizeQuantiles(cbind(l[[30]], l[[31]], l[[32]],l[[33]], l[[34]], l[[35]])) 

mean1=mean(m1)
mean2=mean(m2)
mean3=mean(m3)
mean4=mean(m4)

print(mean1)
print(mean2)
print(mean3)
print(mean4)

#filter regions that are above mean for at least one time point for one 
#histone mark
getThis1 = which((apply(m1,1,function(x) any(x>=mean1)))==TRUE)
getThis2 = which((apply(m2,1,function(x) any(x>=mean2)))==TRUE)
getThis3 = which((apply(m3,1,function(x) any(x>=mean3)))==TRUE)
getThis4 = which((apply(m4,1,function(x) any(x>=mean4)))==TRUE)

m1.filtered=m1[getThis1,]
m2.filtered=m2[getThis2,]
m3.filtered=m3[getThis3,]
m4.filtered=m4[getThis4,]

print(dim(m1.filtered))
print(dim(m2.filtered))
print(dim(m3.filtered))
print(dim(m4.filtered))

g=unique(sort(c(getThis1,getThis2,getThis3,getThis4)))

print(length(g))

filtered1=m1[g,]
filtered2=m2[g,]
filtered3=m3[g,]
filtered4=m4[g,]

print(dim(filtered1))

#write regions and normalized counts for each time point
#it is ordered for D0, D10, D2, D5, D7, D8
#now reorder for each histone modification to D0, D2, D5, D7, D10, D8
#regions 1-11
#histone modification 1: col 12-16 etc
n.reordered = cbind(l[g,seq(1:11)],filtered1[,c(1,3,4,5,2,6)],
                    filtered2[,c(1,3,4,5,2,6)],filtered3[,c(1,3,4,5,2,6)],
                    filtered4[,c(1,3,4,5,2,6)])
print(dim(n.reordered))

write.table(n.reordered, sep = "\t", file = paste0("allCountsNorm.txt"), 
            quote = F, row.names = F, col.names = F)
write.table(l[g,seq(1:11)], sep = "\t", file = paste0("regions_filtered.bed"), 
            quote = F, row.names = F, col.names = F)

f1=m1[g,]+1
f2=m2[g,]+1
f3=m3[g,]+1
f4=m4[g,]+1

#compute FC
#ordering of time points: D0 D2 D5 D7 D10 D8
#ordering of columns: col1 col3 col4 col5 col2 col6
#FC of time points D2/D0 D5/D2 D7/D5 D10/D7 D8/D5
#FC columns: col3/col1 col4/col3 col5/col4 col2/col5 col6/col4

fc1=cbind(log2(f1[,3] / f1[,1]), log2(f1[,4] / f1[,3]), log2(f1[,5] / f1[,4]), 
          log2(f1[,2] / f1[,5], log2(f1[,6] / f1[,4]))
fc2=cbind(log2(f2[,3] / f2[,1]), log2(f2[,4] / f2[,3]), log2(f2[,5] / f2[,4]), 
          log2(f2[,2] / f2[,5], log2(f2[,6] / f2[,4]))
fc3=cbind(log2(f3[,3] / f3[,1]), log2(f3[,4] / f3[,3]), log2(f3[,5] / f3[,4]), 
          log2(f3[,2] / f3[,5], log2(f3[,6] / f3[,4]))
fc4=cbind(log2(f4[,3] / f4[,1]), log2(f4[,4] / f4[,3]), log2(f4[,5] / f4[,4]), 
          log2(f4[,2] / f4[,5], log2(f5[,6] / f5[,4]))

fc.reordered = cbind(l[g,seq(1:11)],fc1*100,fc2*100,fc3*100,fc4*100)

print(dim(fc.reordered))
write.table(fc.reordered, sep = "\t", file = paste0("allFold.txt"), quote = F, 
            row.names = F, col.names = F)

write.table(fc.reordered[,-(1:11)], sep = "\t", 
            file = paste0("allFold_data.txt"), quote = F, row.names = F, 
            col.names = F)

#q()


#histogram of fold changes
#around 0: stays the same
#right side: increasing
#left side: decreasing
#plot distribution of log-FC, should be normal/Gaussian


pdf(paste0("logFC_H3K27ac", name, ".pdf"), width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
hist(fc1[,1], breaks = 100, main=paste0("D2/D0"), xlab=("log2FC"))
hist(fc1[,2], breaks = 100, main=paste0("D5/D2"), xlab=("log2FC"))
hist(fc1[,3], breaks = 100, main=paste0("D7/D5"), xlab=("log2FC"))
hist(fc1[,4], breaks = 100, main=paste0("D10/D7"), xlab=("log2FC"))
hist(fc1[,5], breaks = 100, main=paste0("D8/D5"), xlab=("log2FC"))
mtext("H3K27ac", outer=TRUE, cex=1.5)
dev.off()

pdf(paste0("logFC_H3K27me3", name, ".pdf"), width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
hist(fc2[,1], breaks = 100, main=paste0("D2/D0"), xlab=("log2FC"))
hist(fc2[,2], breaks = 100, main=paste0("D5/D2"), xlab=("log2FC"))
hist(fc2[,3], breaks = 100, main=paste0("D7/D5"), xlab=("log2FC"))
hist(fc2[,4], breaks = 100, main=paste0("D10/D7"), xlab=("log2FC"))
hist(fc2[,5], breaks = 100, main=paste0("D8/D5"), xlab=("log2FC"))
mtext("H3K27me3", outer=TRUE, cex=1.5)
dev.off()

pdf(paste0("logFC_H3K4me1", name, ".pdf"), width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
hist(fc3[,1], breaks = 100, main=paste0("D2/D0"), xlab=("log2FC"))
hist(fc3[,2], breaks = 100, main=paste0("D5/D2"), xlab=("log2FC"))
hist(fc3[,3], breaks = 100, main=paste0("D7/D5"), xlab=("log2FC"))
hist(fc3[,4], breaks = 100, main=paste0("D10/D7"), xlab=("log2FC"))
hist(fc3[,5], breaks = 100, main=paste0("D8/D5"), xlab=("log2FC"))
mtext("H3K4me1", outer=TRUE, cex=1.5)
dev.off()

pdf(paste0("logFC_H3K4me3", name, ".pdf"), width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
hist(fc4[,1], breaks = 100, main=paste0("D2/D0"), xlab=("log2FC"))
hist(fc4[,2], breaks = 100, main=paste0("D5/D2"), xlab=("log2FC"))
hist(fc4[,3], breaks = 100, main=paste0("D7/D5"), xlab=("log2FC"))
hist(fc4[,4], breaks = 100, main=paste0("D10/D7"), xlab=("log2FC"))
hist(fc4[,5], breaks = 100, main=paste0("D8/D5"), xlab=("log2FC"))
mtext("H3K4me3", outer=TRUE, cex=1.5)
dev.off()


#plot distribution of counts

pdf(paste0("Hist_H3K27ac", name, ".pdf"), width=8, height=6)
hist(m1, breaks = 100, main=paste0("H3K27ac"), xlab=("normalized counts"))
lines(rep(mean1, 100), seq(0,200000,length.out=100), col = "red", lty = 1, 
      lwd=2)
dev.off()

pdf(paste0("Hist_H3K27me3", name, ".pdf"), width=8, height=6)
hist(m2, breaks = 100, main=paste0("H3K27me3"), xlab=("normalized counts"))
lines(rep(mean2, 100), seq(0,200000,length.out=100), col = "red", lty = 1, 
      lwd=2)
dev.off()

pdf(paste0("Hist_H3K4me1", name, ".pdf"), width=8, height=6)
hist(m3, breaks = 100, main=paste0("H3K4me1"), xlab=("normalized counts"))
lines(rep(mean3, 100), seq(0,200000,length.out=100), col = "red", lty = 1, 
      lwd=2)
dev.off()

pdf(paste0("Hist_H3K4me3", name, ".pdf"), width=8, height=6)
hist(m4, breaks = 100, main=paste0("H3K4me3"), xlab=("normalized counts"))
lines(rep(mean4, 100), seq(0,200000,length.out=100), col = "red", lty = 1,
      lwd=2)
dev.off()



q()

