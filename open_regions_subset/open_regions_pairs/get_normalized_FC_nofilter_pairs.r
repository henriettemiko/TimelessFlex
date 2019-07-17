##########
#name:          get_normalized_FC_nofilter_pairs.r
#description:   compute normalized fold changes for pairs without filtering
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 16, 2019
##########

library(limma)

args = commandArgs(trailingOnly=TRUE)

curdir <- getwd()
setwd(curdir)

l = read.table(args[1], sep="\t", col.names=1:55)
name <- args[2]

print(dim(l))

#select columns for histone mark (columns = time points)
#regions have 11 columns at beginning of file, 4 time points
#HiC info is col 1, region left is 2-12, counts left 13-28, 
#region right 29-39, counts right 40-55
m1prom = normalizeQuantiles(cbind(l[[13]], l[[14]], l[[15]],l[[16]])) 
m2prom = normalizeQuantiles(cbind(l[[17]], l[[18]], l[[19]],l[[20]])) 
m3prom = normalizeQuantiles(cbind(l[[21]], l[[22]], l[[23]],l[[24]])) 
m4prom = normalizeQuantiles(cbind(l[[25]], l[[26]], l[[27]],l[[28]])) 

m1enh = normalizeQuantiles(cbind(l[[40]], l[[41]], l[[42]],l[[43]])) 
m2enh = normalizeQuantiles(cbind(l[[44]], l[[45]], l[[46]],l[[47]])) 
m3enh = normalizeQuantiles(cbind(l[[48]], l[[49]], l[[50]],l[[51]])) 
m4enh = normalizeQuantiles(cbind(l[[52]], l[[53]], l[[54]],l[[55]])) 

print(head(m1prom))
print(head(m2prom))
print(head(m3prom))
print(head(m4prom))

print(head(m1enh))
print(head(m2enh))
print(head(m3enh))
print(head(m4enh))

mean1prom=mean(m1prom)
mean2prom=mean(m2prom)
mean3prom=mean(m3prom)
mean4prom=mean(m4prom)

mean1enh=mean(m1enh)
mean2enh=mean(m2enh)
mean3enh=mean(m3enh)
mean4enh=mean(m4enh)

mean1=mean(cbind(m1prom,m1enh))
mean2=mean(cbind(m2prom,m2enh))
mean3=mean(cbind(m3prom,m3enh))
mean4=mean(cbind(m4prom,m4enh))


#filter regions that are above mean for at least one time point for one 
#histone mark
#regions not passing are considered noise
good1prom=which((apply(m1prom, 1, function(x) any( x >= mean1prom))) == TRUE)
good2prom=which((apply(m2prom, 1, function(x) any( x >= mean2prom))) == TRUE)
good3prom=which((apply(m3prom, 1, function(x) any( x >= mean3prom))) == TRUE)
good4prom=which((apply(m4prom, 1, function(x) any( x >= mean4prom))) == TRUE)

good1enh= which((apply(m1enh, 1, function(x) any( x >= mean1enh))) == TRUE)
good2enh= which((apply(m2enh, 1, function(x) any( x >= mean2enh))) == TRUE)
good3enh= which((apply(m3enh, 1, function(x) any( x >= mean3enh))) == TRUE)
good4enh= which((apply(m4enh, 1, function(x) any( x >= mean4enh))) == TRUE)

gprom = unique(sort(c(good1prom,good2prom,good3prom,good4prom)))
genh = unique(sort(c(good1enh,good2enh,good3enh,good4enh)))

gboth=sort(c(gprom,genh))
idx = duplicated(gboth)
#g has indices for rows of pairs where prom AND enh have one value above 
#mean for any histone mark at any time
g = gboth[idx]

gnoise = gboth[-idx]

print(length(g))
print(length(gnoise))


m1promfiltered = m1prom[g,]
m2promfiltered = m2prom[g,]
m3promfiltered = m3prom[g,]
m4promfiltered = m4prom[g,]
m1enhfiltered = m1enh[g,]
m2enhfiltered = m2enh[g,]
m3enhfiltered = m3enh[g,]
m4enhfiltered = m4enh[g,]

m1promnoise = m1prom[-g,]
m2promnoise = m2prom[-g,]
m3promnoise = m3prom[-g,]
m4promnoise = m4prom[-g,]
m1enhnoise = m1enh[-g,]
m2enhnoise = m2enh[-g,]
m3enhnoise = m3enh[-g,]
m4enhnoise = m4enh[-g,]


#write regions and normalized counts for each time point
#it is ordered for D0, D10, D2, D5
#now reorder for each histone modification to D0, D2, D5, D10

#HiC info is col 1, region left is 2-12, counts left 13-28, 
#region right 29-39, counts right 40-55

n.reordered.filtered = cbind(l[g,1:12],m1promfiltered[,c(1,3,4,2)],
                             m2promfiltered[,c(1,3,4,2)],
                             m3promfiltered[,c(1,3,4,2)],
                             m4promfiltered[,c(1,3,4,2)],
                             l[g,29:39],m1enhfiltered[,c(1,3,4,2)],
                             m2enhfiltered[,c(1,3,4,2)],
                             m3enhfiltered[,c(1,3,4,2)],
                             m4enhfiltered[,c(1,3,4,2)])

n.reordered.noise = cbind(l[-g,1:12],m1promnoise[,c(1,3,4,2)],
                          m2promnoise[,c(1,3,4,2)],m3promnoise[,c(1,3,4,2)],
                          m4promnoise[,c(1,3,4,2)],l[-g,29:39],
                          m1enhnoise[,c(1,3,4,2)],m2enhnoise[,c(1,3,4,2)],
                          m3enhnoise[,c(1,3,4,2)],m4enhnoise[,c(1,3,4,2)])


write.table(rbind(n.reordered.filtered, n.reordered.noise), 
            sep = "\t", file = paste0("allCountsNorm_", name, ".txt"), 
            quote = F, row.names = F, col.names = F)

write.table(cbind(rbind(l[g,1:12],l[-g,1:12]),
                  rbind(l[g,29:39],l[-g,29:39])), sep = "\t", 
            file = paste0("regions_", name, ".bed"), 
            quote = F, row.names = F, col.names = F)

write.table(n.reordered.filtered, sep = "\t",
            file = paste0("allCountsNorm_filtered_", name, ".txt"), 
            quote = F, row.names = F, col.names = F)

write.table(cbind(l[g,1:12],l[g,29:39]), 
            sep = "\t", file = paste0("regions_filtered_", name, ".bed"), 
            quote = F, row.names = F, col.names = F)

write.table(n.reordered.noise, sep = "\t",
            file = paste0("allCountsNorm_noise_", name, ".txt"), 
            quote = F, row.names = F, col.names = F)

write.table(cbind(l[-g,1:12],l[-g,29:39]), sep = "\t", 
            file = paste0("regions_noise_", name, ".bed"), 
            quote = F, row.names = F, col.names = F)


f1prom=m1promfiltered+1
f2prom=m2promfiltered+1
f3prom=m3promfiltered+1
f4prom=m4promfiltered+1

f1enh=m1enhfiltered+1
f2enh=m2enhfiltered+1
f3enh=m3enhfiltered+1
f4enh=m4enhfiltered+1

f1promnoise=m1promnoise+1
f2promnoise=m2promnoise+1
f3promnoise=m3promnoise+1
f4promnoise=m4promnoise+1

f1enhnoise=m1enhnoise+1
f2enhnoise=m2enhnoise+1
f3enhnoise=m3enhnoise+1
f4enhnoise=m4enhnoise+1


#compute FC
#ordering of time points: D0 D2 D5 D10
#ordering of columns: col1 col3 col4 col2
#FC of time points D2/D0 D5/D2 D10/D5
#FC columns: col3/col1 col4/col3 col2/col4

fc1prom=cbind(log2(f1prom[,3] / f1prom[,1]), log2(f1prom[,4] / f1prom[,3]), 
              log2(f1prom[,2] / f1prom[,4]))
fc2prom=cbind(log2(f2prom[,3] / f2prom[,1]), log2(f2prom[,4] / f2prom[,3]), 
              log2(f2prom[,2] / f2prom[,4]))
fc3prom=cbind(log2(f3prom[,3] / f3prom[,1]), log2(f3prom[,4] / f3prom[,3]), 
              log2(f3prom[,2] / f3prom[,4]))
fc4prom=cbind(log2(f4prom[,3] / f4prom[,1]), log2(f4prom[,4] / f4prom[,3]), 
              log2(f4prom[,2] / f4prom[,4]))


fc1enh=cbind(log2(f1enh[,3] / f1enh[,1]), log2(f1enh[,4] / f1enh[,3]), 
             log2(f1enh[,2] / f1enh[,4]))
fc2enh=cbind(log2(f2enh[,3] / f2enh[,1]), log2(f2enh[,4] / f2enh[,3]), 
             log2(f2enh[,2] / f2enh[,4]))
fc3enh=cbind(log2(f3enh[,3] / f3enh[,1]), log2(f3enh[,4] / f3enh[,3]), 
             log2(f3enh[,2] / f3enh[,4]))
fc4enh=cbind(log2(f4enh[,3] / f4enh[,1]), log2(f4enh[,4] / f4enh[,3]), 
             log2(f4enh[,2] / f4enh[,4]))


fc.reordered.filtered = cbind(l[g,1:12],fc1prom*100,fc2prom*100,fc3prom*100,
                              fc4prom*100,l[g,29:39],fc1enh*100,fc2enh*100,
                              fc3enh*100,fc4enh*100)

fc.reordered.filtered.data = cbind(fc1prom*100,fc2prom*100,fc3prom*100,
                                   fc4prom*100,fc1enh*100,fc2enh*100,
                                   fc3enh*100,fc4enh*100)


print(dim(fc.reordered.filtered))

fc1promnoise=cbind(log2(f1promnoise[,3] / f1promnoise[,1]), 
                   log2(f1promnoise[,4] / f1promnoise[,3]), 
                   log2(f1promnoise[,2] / f1promnoise[,4]))
fc2promnoise=cbind(log2(f2promnoise[,3] / f2promnoise[,1]), 
                   log2(f2promnoise[,4] / f2promnoise[,3]), 
                   log2(f2promnoise[,2] / f2promnoise[,4]))
fc3promnoise=cbind(log2(f3promnoise[,3] / f3promnoise[,1]), 
                   log2(f3promnoise[,4] / f3promnoise[,3]), 
                   log2(f3promnoise[,2] / f3promnoise[,4]))
fc4promnoise=cbind(log2(f4promnoise[,3] / f4promnoise[,1]), 
                   log2(f4promnoise[,4] / f4promnoise[,3]), 
                   log2(f4promnoise[,2] / f4promnoise[,4]))

fc1enhnoise=cbind(log2(f1enhnoise[,3] / f1enhnoise[,1]), 
                  log2(f1enhnoise[,4] / f1enhnoise[,3]), 
                  log2(f1enhnoise[,2] / f1enhnoise[,4]))
fc2enhnoise=cbind(log2(f2enhnoise[,3] / f2enhnoise[,1]), 
                  log2(f2enhnoise[,4] / f2enhnoise[,3]), 
                  log2(f2enhnoise[,2] / f2enhnoise[,4]))
fc3enhnoise=cbind(log2(f3enhnoise[,3] / f3enhnoise[,1]), 
                  log2(f3enhnoise[,4] / f3enhnoise[,3]), 
                  log2(f3enhnoise[,2] / f3enhnoise[,4]))
fc4enhnoise=cbind(log2(f4enhnoise[,3] / f4enhnoise[,1]), 
                  log2(f4enhnoise[,4] / f4enhnoise[,3]), 
                  log2(f4enhnoise[,2] / f4enhnoise[,4]))

fc.reordered.noise = cbind(l[-g,1:12],fc1promnoise*100,
                           fc2promnoise*100,fc3promnoise*100,
                           fc4promnoise*100,l[-g,29:39],
                           fc1enhnoise*100,fc2enhnoise*100,
                           fc3enhnoise*100,fc4enhnoise*100)

fc.reordered.noise.data = cbind(fc1promnoise*100,fc2promnoise*100,
                                fc3promnoise*100,fc4promnoise*100,
                                fc1enhnoise*100, fc2enhnoise*100,
                                fc3enhnoise*100,fc4enhnoise*100)


print(dim(fc.reordered.noise))


write.table(rbind(fc.reordered.filtered, fc.reordered.noise), 
            sep = "\t", file = paste0("allFold_", name, ".txt"), 
            quote = F, row.names = F, col.names = F)

write.table(rbind(fc.reordered.filtered.data, fc.reordered.noise.data), 
            sep = "\t", file = paste0("allFold_data_", name, ".txt"), 
            quote = F, row.names = F, col.names = F)


write.table(fc.reordered.filtered, sep = "\t",
            file = paste0("allFold_filtered_", name, ".txt"), 
            quote = F, row.names = F, col.names = F)

write.table(fc.reordered.noise, sep = "\t",
            file = paste0("allFold_noise_", name, ".txt"), 
            quote = F, row.names = F, col.names = F)







#q()




#histogram of fold changes
#around 0: stays the same
#right side: increasing
#left side: decreasing
#plot distribution of log-FC, should be normal/Gaussian


pdf(paste0("logFC_H3K27ac_", name, ".pdf"), width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
hist(c(fc1prom[,1],fc1enh[,1],fc1promnoise[,1],fc1enhnoise[,1]), 
     breaks = 100, main=paste0("D2/D0"), 
     xlab=("log2FC"))
hist(c(fc1prom[,2],fc1enh[,2],fc1promnoise[,2],fc1enhnoise[,2]), 
     breaks = 100, main=paste0("D5/D2"), 
     xlab=("log2FC"))
hist(c(fc1prom[,3],fc1enh[,3],fc1promnoise[,3],fc1enhnoise[,3]), 
     breaks = 100, main=paste0("D10/D5"), 
     xlab=("log2FC"))
mtext("H3K27ac", outer=TRUE, cex=1.5)
dev.off()

pdf(paste0("logFC_H3K27me3_", name, ".pdf"), width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
hist(c(fc2prom[,1],fc2enh[,1],fc2promnoise[,1],fc2enhnoise[,1]), 
     breaks = 100, main=paste0("D2/D0"), 
     xlab=("log2FC"))
hist(c(fc2prom[,2],fc2enh[,2],fc2promnoise[,2],fc2enhnoise[,2]), 
     breaks = 100, main=paste0("D5/D2"), 
     xlab=("log2FC"))
hist(c(fc2prom[,3],fc2enh[,3],fc2promnoise[,3],fc2enhnoise[,3]), 
     breaks = 100, main=paste0("D10/D5"), 
     xlab=("log2FC"))
mtext("H3K27me3", outer=TRUE, cex=1.5)
dev.off()

pdf(paste0("logFC_H3K4me1_", name, ".pdf"), width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
hist(c(fc3prom[,1],fc3enh[,1],fc3promnoise[,1],fc3enhnoise[,1]), 
     breaks = 100, main=paste0("D2/D0"), 
     xlab=("log2FC"))
hist(c(fc3prom[,2],fc3enh[,2],fc3promnoise[,2],fc3enhnoise[,2]), 
     breaks = 100, main=paste0("D5/D2"), 
     xlab=("log2FC"))
hist(c(fc3prom[,3],fc3enh[,3],fc3promnoise[,3],fc3enhnoise[,3]), 
     breaks = 100, main=paste0("D10/D5"), 
     xlab=("log2FC"))
mtext("H3K4me1", outer=TRUE, cex=1.5)
dev.off()

pdf(paste0("logFC_H3K4me3_", name, ".pdf"), width=8, height=6)
par(mfrow=c(2,2),oma=c(0,0,2,0))
hist(c(fc4prom[,1],fc4enh[,1],fc4promnoise[,1],fc4enhnoise[,1]), 
     breaks = 100, main=paste0("D2/D0"), 
     xlab=("log2FC"))
hist(c(fc4prom[,2],fc4enh[,2],fc4promnoise[,2],fc4enhnoise[,2]), 
     breaks = 100, main=paste0("D5/D2"), 
     xlab=("log2FC"))
hist(c(fc4prom[,3],fc4enh[,3],fc4promnoise[,3],fc4enhnoise[,3]), 
     breaks = 100, main=paste0("D10/D5"), 
     xlab=("log2FC"))
mtext("H3K4me3", outer=TRUE, cex=1.5)
dev.off()


#plot distribution of counts

if (name=="init_prom-enh"){
    left="promoter"
    right="enhancer"
} else if (name=="init_prom-prom"){
    left="promoter"
    right="promoter"
} else if (name=="init_enh-enh"){
    left="enhancer"
    right="enhancer"
}

pdf(paste0("Hist_H3K27ac_", name, ".pdf"), width=8, height=6)
par(mfrow=c(1,2),oma=c(0,0,2,0))
hist(m1prom, breaks = 100, main=paste0(left), xlab=("normalized counts"))
lines(rep(mean1prom, 100), seq(0,200000,length.out=100), col = "red", 
      lty = 1, lwd=2)
hist(m1enh, breaks = 100, main=paste0(right), xlab=("normalized counts"))
lines(rep(mean1enh, 100), seq(0,200000,length.out=100), col = "red", 
      lty = 1, lwd=2)
mtext("H3K27ac", outer=TRUE, cex=1.5)
dev.off()

pdf(paste0("Hist_H3K27me3_", name, ".pdf"), width=8, height=6)
par(mfrow=c(1,2),oma=c(0,0,2,0))
hist(m2prom, breaks = 100, main=paste0(left), xlab=("normalized counts"))
lines(rep(mean2prom, 100), seq(0,200000,length.out=100), col = "red", 
       lty = 1, lwd=2)
hist(m2enh, breaks = 100, main=paste0(right), xlab=("normalized counts"))
lines(rep(mean2enh, 100), seq(0,200000,length.out=100), col = "red", 
      lty = 1, lwd=2)
mtext("H3K27me3", outer=TRUE, cex=1.5)
dev.off()

pdf(paste0("Hist_H3K4me1_", name, ".pdf"), width=8, height=6)
par(mfrow=c(1,2),oma=c(0,0,2,0))
hist(m3prom, breaks = 100, main=paste0(left), xlab=("normalized counts"))
lines(rep(mean3prom, 100), seq(0,200000,length.out=100), col = "red", 
      lty = 1, lwd=2)
hist(m3enh, breaks = 100, main=paste0(right), xlab=("normalized counts"))
lines(rep(mean3enh, 100), seq(0,200000,length.out=100), col = "red", 
      lty = 1, lwd=2)
mtext("H3K4me1", outer=TRUE, cex=1.5)
dev.off()

pdf(paste0("Hist_H3K4me3_", name, ".pdf"), width=8, height=6)
par(mfrow=c(1,2),oma=c(0,0,2,0))
hist(m4prom, breaks = 100, main=paste0(left), xlab=("normalized counts"))
lines(rep(mean4prom, 100), seq(0,200000,length.out=100), col = "red", 
      lty = 1, lwd=2)
hist(m4enh, breaks = 100, main=paste0(right), xlab=("normalized counts"))
lines(rep(mean4enh, 100), seq(0,200000,length.out=100), col = "red", 
      lty = 1, lwd=2)
mtext("H3K4me3", outer=TRUE, cex=1.5)
dev.off()


q()

