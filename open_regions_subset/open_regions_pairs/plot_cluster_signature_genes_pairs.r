
##########
#name:          plot_clusters_signature_genes_pairs.r
#description:   plots clusters and signature genes for pairs
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 29, 2019
##########


library(psych)

setwd(getwd())

args = commandArgs(trailingOnly=TRUE)

numClusters = as.numeric(args[1])
numMarks = as.numeric(args[2])
numTimePoints = as.numeric(args[3])
timeless.dir=args[4]
signalgenerator.dir=args[5]
numcluster.dir=args[6]
type=args[7]

score <- function(x){
    x = scale(x, scale = FALSE)
    x  = ((x-min(x))/(max(x)-min(x))) * 100
    return(x)
}


filename=""
if (type=="multi_prom-enh" || type=="multi_prom-prom" || 
    type=="multi_enh-enh"){
    filename="_afterEM"
} 
print(filename)

filepath=""
#if (type=="init_prom-enh" || type=="init_prom-prom" || 
#    type=="init_enh-enh"){
#    filepath="/../"
#} 
print(filepath)


l = read.table(paste0(timeless.dir,filepath,"/classes-", numClusters, 
                      filename, ".txt"))
k = read.table(paste0(signalgenerator.dir,"/allCountsNorm.txt")) 


#score separately
#PROMOTER timepoints, ENHANCER timepoints
# 4 time points prom, 4 timepoints enh
k27acprom = cbind(k[[13]], k[[14]], k[[15]], k[[16]])
k27me3prom = cbind(k[[17]], k[[18]], k[[19]], k[[20]])
k4me1prom = cbind(k[[21]], k[[22]], k[[23]], k[[24]])
k4me3prom = cbind(k[[25]], k[[26]], k[[27]], k[[28]])

k27acenh = cbind(k[[40]], k[[41]], k[[42]], k[[43]])
k27me3enh = cbind(k[[44]], k[[45]], k[[46]], k[[47]])
k4me1enh = cbind(k[[48]], k[[49]], k[[50]], k[[51]])
k4me3enh = cbind(k[[52]], k[[53]], k[[54]], k[[55]])

k27me3prom = matrix(score(as.vector(k27me3prom)), ncol = 4, byrow = FALSE)
k27acprom = matrix(score(as.vector(k27acprom)), ncol = 4, byrow = FALSE)
k4me1prom = matrix(score(as.vector(k4me1prom)), ncol = 4, byrow = FALSE)
k4me3prom = matrix(score(as.vector(k4me3prom)), ncol = 4, byrow = FALSE)

k27me3enh = matrix(score(as.vector(k27me3enh)), ncol = 4, byrow = FALSE)
k27acenh = matrix(score(as.vector(k27acenh)), ncol = 4, byrow = FALSE)
k4me1enh = matrix(score(as.vector(k4me1enh)), ncol = 4, byrow = FALSE)
k4me3enh = matrix(score(as.vector(k4me3enh)), ncol = 4, byrow = FALSE)

k27ac = cbind(k27acprom, k27acenh)
k27me3 = cbind(k27me3prom, k27me3enh)
k4me1 = cbind(k4me1prom, k4me1enh)
k4me3 = cbind(k4me3prom, k4me3enh)

print(colMeans(k27ac))
print(colMeans(k27me3))
print(colMeans(k4me1))
print(colMeans(k4me3))

#read in information about signature genes in clusters
DE = read.table(paste0(numcluster.dir,"/DE_genes_clusters.txt"))
GT = read.table(paste0(numcluster.dir,"/GT_genes_clusters.txt"))
PE = read.table(paste0(numcluster.dir,"/PE_genes_clusters.txt"))

print(DE)
print(GT)
print(PE)

#all signature genes from list that were unambiguously assigned in this 
#clustering
DE.all = read.table(paste0(numcluster.dir,"/DE_genes.txt"))
print(nrow(DE.all))

GT.all = read.table(paste0(numcluster.dir,"/GT_genes.txt"))
print(nrow(GT.all))

PE.all = read.table(paste0(numcluster.dir,"/PE_genes.txt"))
print(nrow(PE.all))


cbbPalette <- c("gray", "#009E73", "#D50F25", "black")
#me3 gray, me1 black, 27ac green 009E73, 27me3 red
#plot from dark to light colors: me1, k27me3, k27ac, me3

pdf(paste0(numcluster.dir,"/clusters_", numClusters, "_signature_genes.pdf"), 
    height=6, width=10) 

par(mar=c(7.1, 4.1, 4.1, 7.1), xpd=TRUE)

par(mfrow = c(1,2))

for (i in 1:numClusters) {

    cur.items = which(l[[1]] == i)
    stuff.k27ac.prom <- cbind(k27ac[cur.items,1], k27ac[cur.items,2], 
                              k27ac[cur.items,3], k27ac[cur.items,4])
    stuff.k27me3.prom <- cbind(k27me3[cur.items,1], k27me3[cur.items,2], 
                               k27me3[cur.items,3], k27me3[cur.items,4])
    stuff.k4me1.prom <- cbind(k4me1[cur.items,1], k4me1[cur.items,2], 
                              k4me1[cur.items,3], k4me1[cur.items,4])
    stuff.k4me3.prom <- cbind(k4me3[cur.items,1], k4me3[cur.items,2], 
                              k4me3[cur.items,3], k4me3[cur.items,4])

    colnames(stuff.k4me1.prom) <- c("D0", "D2", "D5", "D10")



    stuff.k27ac.enh <- cbind(k27ac[cur.items,5], k27ac[cur.items,6], 
                             k27ac[cur.items,7], k27ac[cur.items,8])
    stuff.k27me3.enh <- cbind(k27me3[cur.items,5], k27me3[cur.items,6], 
                              k27me3[cur.items,7], k27me3[cur.items,8])
    stuff.k4me1.enh <- cbind(k4me1[cur.items,5], k4me1[cur.items,6], 
                             k4me1[cur.items,7], k4me1[cur.items,8])
    stuff.k4me3.enh <- cbind(k4me3[cur.items,5], k4me3[cur.items,6], 
                             k4me3[cur.items,7], k4me3[cur.items,8])

    colnames(stuff.k4me1.enh) <- c("D0", "D2", "D5", "D10")


    if (type=="init_prom-enh" || type=="multi_prom-enh"){
        left="PROMOTER"
        right="ENHANCER"
    } else if (type=="init_prom-prom" || type=="multi_prom-prom"){
        left="PROMOTER"
        right="PROMOTER"
    } else if (type=="init_enh-enh" || type=="multi_enh-enh"){
        left="ENHANCER"
        right="ENHANCER"
    }




    #plot promoter

    error.bars(stuff.k4me1.prom, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[4], labels = colnames(stuff.k4me1.prom), 
               ylim = c(0,65), 
               ylab = "normalized average histone mark signal", 
               xlab="time", col = cbbPalette[4], 
               main = paste0("Cluster ", i, " ", left))
    lines(colMeans((cbind(k4me1[cur.items,1], k4me1[cur.items,2], 
                          k4me1[cur.items,3], k4me1[cur.items,4]))), 
          col = cbbPalette[4], ylim = c(0,65), lwd = 8, type="o", 
          pch=16)

    error.bars(stuff.k4me3.prom, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[1], ylim = c(0,65), 
               col = cbbPalette[1], add = TRUE)
    lines(colMeans((cbind(k4me3[cur.items,1], k4me3[cur.items,2], 
                          k4me3[cur.items,3], k4me3[cur.items,4]))), 
          col = cbbPalette[1], ylim = c(0,65), lwd = 8, type="o", 
          pch=16)

    error.bars(stuff.k27ac.prom, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[2], ylim = c(0,65), 
               col = cbbPalette[2], add = TRUE)
    lines(colMeans((cbind(k27ac[cur.items,1], k27ac[cur.items,2], 
                          k27ac[cur.items,3], k27ac[cur.items,4]))), 
          col = cbbPalette[2], ylim = c(0,65), lwd = 8, type="o", 
          pch=16)

    error.bars(stuff.k27me3.prom, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[3], ylim = c(0,65), 
               col = cbbPalette[3], add = TRUE)
    lines(colMeans((cbind(k27me3[cur.items,1], k27me3[cur.items,2], 
                          k27me3[cur.items,3], k27me3[cur.items,4]))), 
          col = cbbPalette[3], ylim = c(0,65), lwd = 8, type="o", 
          pch=16)


    DE.num=DE[DE$V2==i,]$V1
    if (length(DE.num)==0) {
        DE.num = 0
    }

    GT.num=GT[GT$V2==i,]$V1
    if (length(GT.num)==0) {
        GT.num = 0
    }

    PE.num=PE[PE$V2==i,]$V1
    if (length(PE.num)==0) {
        PE.num = 0
    }

    print(DE.num)
    print(GT.num)
    print(PE.num)

    legend("topright", inset=c(-0.5,0), 
           legend=c("H3K27ac","H3K27me3","H3K4me1","H3K4me3"), 
           col=c("#009E73", "#D50F25", "black", "gray"), pch=15) 
    mtext(side=1, text=paste0("signature genes:\nD2: ", DE.num, 
                              " / ", nrow(DE.all), ", D5: ", GT.num, 
                              " / ", nrow(GT.all), 
                              ", D10: ", PE.num, 
                              " / ", nrow(PE.all)), line=6)



    #plot enhancer

    error.bars(stuff.k4me1.enh, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[4], labels = colnames(stuff.k4me1.enh), 
               ylim = c(0,65), 
               ylab = "normalized average histone mark signal", 
               xlab="time", col = cbbPalette[4], 
               main = paste0(right, " (", length(cur.items), " pairs)"))
    lines(colMeans((cbind(k4me1[cur.items,5], k4me1[cur.items,6], 
                          k4me1[cur.items,7], k4me1[cur.items,8]))), 
          col = cbbPalette[4], ylim = c(0,65), lwd = 8, type="o", 
          pch=16)

    error.bars(stuff.k4me3.enh, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[1], ylim = c(0,65), 
               col = cbbPalette[1], add = TRUE)
    lines(colMeans((cbind(k4me3[cur.items,5], k4me3[cur.items,6], 
                          k4me3[cur.items,7], k4me3[cur.items,8]))), 
          col = cbbPalette[1], ylim = c(0,65), lwd = 8, type="o", 
          pch=16)

    error.bars(stuff.k27ac.enh, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[2], ylim = c(0,65), 
               col = cbbPalette[2], add = TRUE)
    lines(colMeans((cbind(k27ac[cur.items,5], k27ac[cur.items,6], 
                          k27ac[cur.items,7], k27ac[cur.items,8]))), 
          col = cbbPalette[2], ylim = c(0,65), lwd = 8, type="o", 
          pch=16)

    error.bars(stuff.k27me3.enh, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[3], ylim = c(0,65), 
               col = cbbPalette[3], add = TRUE)
    lines(colMeans((cbind(k27me3[cur.items,5], k27me3[cur.items,6], 
                          k27me3[cur.items,7], k27me3[cur.items,8]))), 
          col = cbbPalette[3], ylim = c(0,65), lwd = 8, type="o", 
          pch=16)


}


dev.off()
warnings()


q()

