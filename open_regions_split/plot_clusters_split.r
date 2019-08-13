
##########
#name:          plot_clusters_fullset.r
#description:   plots clusters for fullset
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 10, 2019
##########


library(psych)

setwd(getwd())

args = commandArgs(trailingOnly=TRUE)

numClusters = as.numeric(args[1])
numMarks = as.numeric(args[2])
numTimePoints = as.numeric(args[3])
numcluster.dir = args[4]
timeless.dir = args[5]
signalgenerator.dir = args[6]
type = args[7]

score <- function(x){
    x = scale(x, scale = FALSE)
    x  = ((x-min(x))/(max(x)-min(x))) * 100
    return(x)
}

l = read.table(paste0(numcluster.dir, "/classes-", numClusters, ".txt"))
k = read.table(paste0(signalgenerator.dir, "/allCountsNorm.txt"))

# 6 time points
k27ac = cbind(k[[12]], k[[13]], k[[14]], k[[15]], k[[16]], k[[17]])
k27me3 = cbind(k[[18]], k[[19]], k[[20]], k[[21]], k[[22]], k[[23]])
k4me1 = cbind(k[[24]], k[[25]], k[[26]], k[[27]], k[[28]], k[[29]])
k4me3 = cbind(k[[30]], k[[31]], k[[32]], k[[33]], k[[34]], k[[35]])

k27ac = matrix(score(as.vector(k27ac)), ncol = 6, byrow = FALSE)
k27me3 = matrix(score(as.vector(k27me3)), ncol = 6, byrow = FALSE)
k4me1 = matrix(score(as.vector(k4me1)), ncol = 6, byrow = FALSE)
k4me3 = matrix(score(as.vector(k4me3)), ncol = 6, byrow = FALSE)

print(colMeans(k27ac))
print(colMeans(k27me3))
print(colMeans(k4me1))
print(colMeans(k4me3))

cbbPalette <- c("gray", "#009E73", "#D50F25", "black")
#me3 gray, me1 black, 27ac green 009E73, 27me3 red
#plot from dark to light colors: me1, k27me3, k27ac, me3

pdf(paste0(type, "_clusters_", numClusters, ".pdf"), height=6, width=4)

for (i in 1:numClusters) {

    cur.items = which(l[[1]] == i)
    stuff.k27ac <- cbind(k27ac[cur.items,1], k27ac[cur.items,2], 
                         k27ac[cur.items,3], k27ac[cur.items,4], 
                         k27ac[cur.items,5])

    stuff.k27ac.D8 <- cbind(k27ac[cur.items,1], k27ac[cur.items,2], 
                         k27ac[cur.items,3], k27ac[cur.items,6])

    stuff.k27me3 <- cbind(k27me3[cur.items,1], k27me3[cur.items,2], 
                          k27me3[cur.items,3], k27me3[cur.items,4], 
                          k27me3[cur.items,5])

    stuff.k27me3.D8 <- cbind(k27me3[cur.items,1], k27me3[cur.items,2], 
                          k27me3[cur.items,3], k27me3[cur.items,6])

    stuff.k4me1 <- cbind(k4me1[cur.items,1], k4me1[cur.items,2], 
                         k4me1[cur.items,3], k4me1[cur.items,4], 
                         k4me1[cur.items,5])

    stuff.k4me1.D8 <- cbind(k4me1[cur.items,1], k4me1[cur.items,2], 
                         k4me1[cur.items,3], k4me1[cur.items,6])

    stuff.k4me3 <- cbind(k4me3[cur.items,1], k4me3[cur.items,2], 
                         k4me3[cur.items,3], k4me3[cur.items,4], 
                         k4me3[cur.items,5])

    stuff.k4me3.D8 <- cbind(k4me3[cur.items,1], k4me3[cur.items,2], 
                         k4me3[cur.items,3], k4me3[cur.items,6])


    colnames(stuff.k4me1) <- c("D0", "D2", "D5", "D7", "D10")

    colnames(stuff.k4me1.D8) <- c("D0", "D2", "D5", "D8")

    error.bars(stuff.k4me1, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[4], labels = colnames(stuff.k4me1), 
               ylim = c(0,70), 
               ylab = "normalized average histone mark signal", 
               xlab="", col = cbbPalette[4], 
               main = paste0("Cluster ", i, ": n=", length(cur.items)))

    error.bars(stuff.k4me1.D8, add=T, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[4], labels = colnames(stuff.k4me1.D8), 
               ylim = c(0,70), 
               ylab = "normalized average histone mark signal", 
               xlab="", col = cbbPalette[4], 
               main = paste0("Cluster ", i, ": n=", length(cur.items)))

    lines(colMeans((cbind(k4me1[cur.items,1], k4me1[cur.items,2], 
                          k4me1[cur.items,3], k4me1[cur.items,4], 
                          k4me1[cur.items,5]))), col = cbbPalette[4], 
          ylim = c(0,70), lwd = 8, type="o", pch=16)

    lines(colMeans((cbind(k4me1[cur.items,1], k4me1[cur.items,2], 
                          k4me1[cur.items,3], k4me1[cur.items,6]))), 
          col = cbbPalette[4], ylim = c(0,70), lwd = 8, type="o", pch=16)


    error.bars(stuff.k4me3, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[1], ylim = c(0,70), 
               col = cbbPalette[1], add = TRUE)

    error.bars(stuff.k4me3.D8, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[1], ylim = c(0,70), 
               col = cbbPalette[1], add = TRUE)

    lines(colMeans((cbind(k4me3[cur.items,1], k4me3[cur.items,2], 
                          k4me3[cur.items,3], k4me3[cur.items,4], 
                          k4me3[cur.items,5]))), col = cbbPalette[1], 
          ylim = c(0,70), lwd = 8, type="o", pch=16)

    lines(colMeans((cbind(k4me3[cur.items,1], k4me3[cur.items,2], 
                          k4me3[cur.items,3], k4me3[cur.items,6]))), 
          col = cbbPalette[1], ylim = c(0,70), lwd = 8, type="o", pch=16)


    error.bars(stuff.k27ac, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[2], ylim = c(0,70), 
               col = cbbPalette[2], add = TRUE)

    error.bars(stuff.k27ac.D8, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[2], ylim = c(0,70), 
               col = cbbPalette[2], add = TRUE)

    lines(colMeans((cbind(k27ac[cur.items,1], k27ac[cur.items,2], 
                          k27ac[cur.items,3], k27ac[cur.items,4], 
                          k27ac[cur.items,5]))), col = cbbPalette[2], 
          ylim = c(0,70), lwd = 8, type="o", pch=16)

    lines(colMeans((cbind(k27ac[cur.items,1], k27ac[cur.items,2], 
                          k27ac[cur.items,3], k27ac[cur.items,6]))), 
          col = cbbPalette[2], ylim = c(0,70), lwd = 8, type="o", pch=16)


    error.bars(stuff.k27me3, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[3], ylim = c(0,70), 
               col = cbbPalette[3], add = TRUE)

    error.bars(stuff.k27me3.D8, eyes = FALSE, sd = FALSE, bars = FALSE, 
               arrow.col = cbbPalette[3], ylim = c(0,70), 
               col = cbbPalette[3], add = TRUE)

    lines(colMeans((cbind(k27me3[cur.items,1], k27me3[cur.items,2], 
                          k27me3[cur.items,3], k27me3[cur.items,4], 
                          k27me3[cur.items,5]))), col = cbbPalette[3], 
          ylim = c(0,70), lwd = 8, type="o", pch=16)

    lines(colMeans((cbind(k27me3[cur.items,1], k27me3[cur.items,2], 
                          k27me3[cur.items,3], k27me3[cur.items,6]))), 
          col = cbbPalette[3], ylim = c(0,70), lwd = 8, type="o", pch=16)


    axis(side=1, labels=c("D0", "D2", "D5", "D8"), at = c(1,2,3,4), line=2)
    mtext(side=1, "time", line=4)
  


}


dev.off()
warnings()


q()

