##########
#name:          get_normalized_FC_nofilter.r
#description:   compute normalized fold changes for subset without filtering
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 12, 2019
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
#regions have 11 columns at beginning of file, 4 time points
m1 = normalizeQuantiles(cbind(l[[12]], l[[13]], l[[14]],l[[15]])) 
m2 = normalizeQuantiles(cbind(l[[16]], l[[17]], l[[18]],l[[19]])) 
m3 = normalizeQuantiles(cbind(l[[20]], l[[21]], l[[22]],l[[23]])) 
m4 = normalizeQuantiles(cbind(l[[24]], l[[25]], l[[26]],l[[27]])) 

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
print(length(-g))

filtered1=m1[g,]
filtered2=m2[g,]
filtered3=m3[g,]
filtered4=m4[g,]

#regions not passing are considered noise
noise1=m1[-g,]
noise2=m2[-g,]
noise3=m3[-g,]
noise4=m4[-g,]




#write regions and normalized counts for each time point
#it is ordered for D0, D10, D2, D5
#now reorder for each histone modification to D0, D2, D5, D10
#regions 1-11
#histone modification 1: col 12-15 etc



n.reordered.filtered = cbind(l[g,seq(1:11)],filtered1[,c(1,3,4,2)],
                    filtered2[,c(1,3,4,2)],filtered3[,c(1,3,4,2)],
                    filtered4[,c(1,3,4,2)])

n.reordered.noise = cbind(l[-g,seq(1:11)],noise1[,c(1,3,4,2)],
                    noise2[,c(1,3,4,2)],noise3[,c(1,3,4,2)],
                    noise4[,c(1,3,4,2)])

print(dim(n.reordered.filtered))
print(dim(n.reordered.noise))

write.table(rbind(n.reordered.filtered, n.reordered.noise), 
            sep = "\t", file = paste0("allCountsNorm.txt"), 
            quote = F, row.names = F, col.names = F)

write.table(rbind(l[g,seq(1:11)],l[-g,seq(1:11)]), sep = "\t", 
            file = paste0("regions.bed"), 
            quote = F, row.names = F, col.names = F)

write.table(n.reordered.filtered, sep = "\t",
            file = paste0("allCountsNorm_filtered.txt"), 
            quote = F, row.names = F, col.names = F)

write.table(l[g,seq(1:11)], sep = "\t", file = paste0("regions_filtered.bed"), 
            quote = F, row.names = F, col.names = F)

write.table(n.reordered.noise, sep = "\t",
            file = paste0("allCountsNorm_noise.txt"), 
            quote = F, row.names = F, col.names = F)

write.table(l[-g,seq(1:11)], sep = "\t", file = paste0("regions_noise.bed"), 
            quote = F, row.names = F, col.names = F)



f1=filtered1+1
f2=filtered2+1
f3=filtered3+1
f4=filtered4+1

n1=noise1+1
n2=noise2+1
n3=noise3+1
n4=noise4+1

#compute FC
#ordering of time points: D0 D2 D5 D10
#ordering of columns: col1 col3 col4 col2
#FC of time points D2/D0 D5/D2 D10/D5
#FC columns: col3/col1 col4/col3 col2/col4

fc1=cbind(log2(f1[,3] / f1[,1]), log2(f1[,4] / f1[,3]), log2(f1[,2] / f1[,4]))
fc2=cbind(log2(f2[,3] / f2[,1]), log2(f2[,4] / f2[,3]), log2(f2[,2] / f2[,4]))
fc3=cbind(log2(f3[,3] / f3[,1]), log2(f3[,4] / f3[,3]), log2(f3[,2] / f3[,4]))
fc4=cbind(log2(f4[,3] / f4[,1]), log2(f4[,4] / f4[,3]), log2(f4[,2] / f4[,4]))

fc.reordered = cbind(l[g,seq(1:11)],fc1*100,fc2*100,fc3*100,fc4*100)

print(dim(fc.reordered))


nfc1=cbind(log2(n1[,3] / n1[,1]), log2(n1[,4] / n1[,3]), log2(n1[,2] / n1[,4]))
nfc2=cbind(log2(n2[,3] / n2[,1]), log2(n2[,4] / n2[,3]), log2(n2[,2] / n2[,4]))
nfc3=cbind(log2(n3[,3] / n3[,1]), log2(n3[,4] / n3[,3]), log2(n3[,2] / n3[,4]))
nfc4=cbind(log2(n4[,3] / n4[,1]), log2(n4[,4] / n4[,3]), log2(n4[,2] / n4[,4]))

nfc.reordered = cbind(l[-g,seq(1:11)],nfc1*100,nfc2*100,nfc3*100,nfc4*100)

print(dim(nfc.reordered))


write.table(rbind(fc.reordered, nfc.reordered), 
            sep = "\t", file = paste0("allFold.txt"), 
            quote = F, row.names = F, col.names = F)

write.table(rbind(fc.reordered[,-(1:11)],nfc.reordered[,-(1:11)]), sep = "\t", 
            file = paste0("allFold_data.txt"), quote = F, row.names = F, 
            col.names = F)

write.table(fc.reordered, sep = "\t", file = paste0("allFold_filtered.txt"), 
            quote = F, row.names = F, col.names = F)

write.table(nfc.reordered, sep = "\t", file = paste0("allFold_noise.txt"), 
            quote = F, row.names = F, col.names = F)


#q()


#histogram of fold changes
#around 0: stays the same
#right side: increasing
#left side: decreasing
#plot distribution of log-FC, should be normal/Gaussian


pdf(paste0("logFC_H3K27ac", name, ".pdf"), width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
hist(c(fc1[,1],nfc1[,1]), breaks = 100, main=paste0("D2/D0"), 
     xlab=("log2FC"))
hist(c(fc1[,2],nfc1[,2]), breaks = 100, main=paste0("D5/D2"), 
     xlab=("log2FC"))
hist(c(fc1[,3],nfc1[,3]), breaks = 100, main=paste0("D10/D5"), 
     xlab=("log2FC"))
mtext("H3K27ac", outer=TRUE, cex=1.5)
dev.off()

pdf(paste0("logFC_H3K27me3", name, ".pdf"), width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
hist(c(fc2[,1],nfc2[,1]), breaks = 100, main=paste0("D2/D0"), 
     xlab=("log2FC"))
hist(c(fc2[,2],nfc2[,2]), breaks = 100, main=paste0("D5/D2"), 
     xlab=("log2FC"))
hist(c(fc2[,3],nfc2[,3]), breaks = 100, main=paste0("D10/D5"), 
     xlab=("log2FC"))
mtext("H3K27me3", outer=TRUE, cex=1.5)
dev.off()

pdf(paste0("logFC_H3K4me1", name, ".pdf"), width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
hist(c(fc3[,1],nfc3[,1]), breaks = 100, main=paste0("D2/D0"), 
     xlab=("log2FC"))
hist(c(fc3[,2],nfc3[,2]), breaks = 100, main=paste0("D5/D2"), 
     xlab=("log2FC"))
hist(c(fc3[,3],nfc3[,3]), breaks = 100, main=paste0("D10/D5"), 
     xlab=("log2FC"))
mtext("H3K4me1", outer=TRUE, cex=1.5)
dev.off()

pdf(paste0("logFC_H3K4me3", name, ".pdf"), width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
hist(c(fc4[,1],nfc4[,1]), breaks = 100, main=paste0("D2/D0"), 
     xlab=("log2FC"))
hist(c(fc4[,2],nfc4[,2]), breaks = 100, main=paste0("D5/D2"), 
     xlab=("log2FC"))
hist(c(fc4[,3],nfc4[,3]), breaks = 100, main=paste0("D10/D5"), 
     xlab=("log2FC"))
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

