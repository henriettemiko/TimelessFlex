
##########
#name:          plot_clusters.r
#description:   plots clusters
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          April 4, 2020
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

#regions col1-11
#6 timepoints: CMP, MEP, EryA, GMP, Granu, Mono
k27ac = cbind(k[[12]], k[[13]], k[[14]], k[[15]], k[[16]], k[[17]])
k4me1 = cbind(k[[18]], k[[19]], k[[20]], k[[21]], k[[22]], k[[23]])
k4me2 = cbind(k[[24]], k[[25]], k[[26]], k[[27]], k[[28]], k[[29]])
k4me3 = cbind(k[[30]], k[[31]], k[[32]], k[[33]], k[[34]], k[[35]])

k27ac = matrix(score(as.vector(k27ac)), ncol = 6, byrow = FALSE)
k4me1 = matrix(score(as.vector(k4me1)), ncol = 6, byrow = FALSE)
k4me2 = matrix(score(as.vector(k4me2)), ncol = 6, byrow = FALSE)
k4me3 = matrix(score(as.vector(k4me3)), ncol = 6, byrow = FALSE)

print(range(k27ac))
print(range(k4me1))
print(range(k4me2))
print(range(k4me3))

print(colMeans(k27ac))
print(colMeans(k4me1))
print(colMeans(k4me2))
print(colMeans(k4me3))

cbbPalette <- c("#009E73", "black", "orange", "gray")
#me3 gray, me1 black, 27ac green 009E73, k4me2 orange
#I want to plot from dark to light colors
#k4me1, k27ac, k4me2, k4me3


pdf(paste0(type, "_clusters_", numClusters, ".pdf"), height=6, width=6)

par(mar=c(5.1, 4.1, 4.1, 8.1), oma=c(0,0,0,0), xpd=TRUE)

for (i in 1:numClusters) {


    
  cur.items = which(l[[1]] == i)
  stuff.k27ac.1 <- cbind(k27ac[cur.items,1:3])
  stuff.k4me1.1 <- cbind(k4me1[cur.items,1:3])
  stuff.k4me2.1 <- cbind(k4me2[cur.items,1:3])
  stuff.k4me3.1 <- cbind(k4me3[cur.items,1:3])

  stuff.k27ac.2 <- cbind(k27ac[cur.items,c(1,4:5)])
  stuff.k4me1.2 <- cbind(k4me1[cur.items,c(1,4:5)])
  stuff.k4me2.2 <- cbind(k4me2[cur.items,c(1,4:5)])
  stuff.k4me3.2 <- cbind(k4me3[cur.items,c(1,4:5)])


  stuff.k27ac.3 <- cbind(k27ac[cur.items,c(1,4,6)])
  stuff.k4me1.3 <- cbind(k4me1[cur.items,c(1,4,6)])
  stuff.k4me2.3 <- cbind(k4me2[cur.items,c(1,4,6)])
  stuff.k4me3.3 <- cbind(k4me3[cur.items,c(1,4,6)])

  colnames(stuff.k27ac.1) <- c("CMP","MEP","EryA")
  colnames(stuff.k4me1.1) <- c("CMP","MEP","EryA")
  colnames(stuff.k4me2.1) <- c("CMP","MEP","EryA")
  colnames(stuff.k4me3.1) <- c("CMP","MEP","EryA")

  colnames(stuff.k27ac.2) <- c("CMP","GMP","Granu")
  colnames(stuff.k4me1.2) <- c("CMP","GMP","Granu")
  colnames(stuff.k4me2.2) <- c("CMP","GMP","Granu")
  colnames(stuff.k4me3.2) <- c("CMP","GMP","Granu")

  colnames(stuff.k27ac.2) <- c("CMP","GMP","Mono")
  colnames(stuff.k4me1.2) <- c("CMP","GMP","Mono")
  colnames(stuff.k4me2.2) <- c("CMP","GMP","Mono")
  colnames(stuff.k4me3.2) <- c("CMP","GMP","Mono")

  #solid lines represent CMP-MEP-EryA, dashed lines CMP-GMP-Granu and dotted lines CMP-GMP-Mono differentiation


#K4me1

  error.bars(stuff.k4me1.1, eyes = FALSE, sd = FALSE, bars = FALSE, arrow.col = cbbPalette[2], labels = colnames(stuff.k4me1.1), ylim = c(0,60), ylab = "normalized average histone mark signal", xlab="", col = cbbPalette[2], main = paste0("Cluster ", i, " (", length(cur.items), " regions)"))
   lines(colMeans(stuff.k4me1.1), col = cbbPalette[2], ylim = c(0,60), lwd = 8, type="o", pch=15)

  error.bars(stuff.k4me1.2, eyes = FALSE, sd = FALSE, bars = FALSE, arrow.col = cbbPalette[2], labels = colnames(stuff.k4me1.2), ylim = c(0,60), col = cbbPalette[2], add=T)
   lines(colMeans(stuff.k4me1.2), col = cbbPalette[2], ylim = c(0,60), lwd = 4, type="o", pch=16, lty=2)

  error.bars(stuff.k4me1.3, eyes = FALSE, sd = FALSE, bars = FALSE, arrow.col = cbbPalette[2], labels = colnames(stuff.k4me1.3), ylim = c(0,60), col = cbbPalette[2], add=T)
   lines(colMeans(stuff.k4me1.3), col = cbbPalette[2], ylim = c(0,60), lwd = 10, type="o", pch=17,lty=3)


#K4me2

  error.bars(stuff.k4me2.1, eyes = FALSE, sd = FALSE, bars = FALSE, arrow.col = cbbPalette[3], labels = colnames(stuff.k4me2.1), ylim = c(0,60), ylab = "normalized average histone mark signal", xlab="", col = cbbPalette[3], add=T)
   lines(colMeans(stuff.k4me2.1), col = cbbPalette[3], ylim = c(0,60), lwd = 8, type="o", pch=15)

  error.bars(stuff.k4me2.2, eyes = FALSE, sd = FALSE, bars = FALSE, arrow.col = cbbPalette[3], labels = colnames(stuff.k4me2.2), ylim = c(0,60), col = cbbPalette[3], add=T)
   lines(colMeans(stuff.k4me2.2), col = cbbPalette[3], ylim = c(0,60), lwd = 4, type="o", pch=16, lty=2)

  error.bars(stuff.k4me2.3, eyes = FALSE, sd = FALSE, bars = FALSE, arrow.col = cbbPalette[3], labels = colnames(stuff.k4me2.3), ylim = c(0,60), col = cbbPalette[3], add=T)
   lines(colMeans(stuff.k4me2.3), col = cbbPalette[3], ylim = c(0,60), lwd = 10, type="o", pch=17,lty=3)


#K4me3

  error.bars(stuff.k4me3.1, eyes = FALSE, sd = FALSE, bars = FALSE, arrow.col = cbbPalette[4], labels = colnames(stuff.k4me3.1), ylim = c(0,60), ylab = "normalized average histone mark signal", xlab="", col = cbbPalette[4], add=T)
   lines(colMeans(stuff.k4me3.1), col = cbbPalette[4], ylim = c(0,60), lwd = 8, type="o", pch=15)

  error.bars(stuff.k4me3.2, eyes = FALSE, sd = FALSE, bars = FALSE, arrow.col = cbbPalette[4], labels = colnames(stuff.k4me3.2), ylim = c(0,60), col = cbbPalette[4], add=T)
   lines(colMeans(stuff.k4me3.2), col = cbbPalette[4], ylim = c(0,60), lwd = 4, type="o", pch=16, lty=2)

  error.bars(stuff.k4me3.3, eyes = FALSE, sd = FALSE, bars = FALSE, arrow.col = cbbPalette[4], labels = colnames(stuff.k4me3.3), ylim = c(0,60), col = cbbPalette[4], add=T)
   lines(colMeans(stuff.k4me3.3), col = cbbPalette[4], ylim = c(0,60), lwd = 10, type="o", pch=17,lty=3)


#K27ac

  error.bars(stuff.k27ac.1, eyes = FALSE, sd = FALSE, bars = FALSE, arrow.col = cbbPalette[1], labels = colnames(stuff.k27ac.1), ylim = c(0,60), ylab = "normalized average histone mark signal", xlab="", col = cbbPalette[1], add=T)
   lines(colMeans(stuff.k27ac.1), col = cbbPalette[1], ylim = c(0,60), lwd = 8, type="o", pch=15)

  error.bars(stuff.k27ac.2, eyes = FALSE, sd = FALSE, bars = FALSE, arrow.col = cbbPalette[1], labels = colnames(stuff.k27ac.2), ylim = c(0,60), col = cbbPalette[1], add=T)
   lines(colMeans(stuff.k27ac.2), col = cbbPalette[1], ylim = c(0,60), lwd = 4, type="o", pch=16, lty=2)

  error.bars(stuff.k27ac.3, eyes = FALSE, sd = FALSE, bars = FALSE, arrow.col = cbbPalette[1], labels = colnames(stuff.k27ac.3), ylim = c(0,60), col = cbbPalette[1], add=T)
   lines(colMeans(stuff.k27ac.3), col = cbbPalette[1], ylim = c(0,60), lwd = 10, type="o", pch=17,lty=3)


axis(side=1, labels=c("","GMP","Granu/Mono"), at = c(1,2,3), line=2)
mtext(side=1, "time", line=4)

legend("topright", inset=c(-0.45,0),
           legend=c("H3K27ac", "H3K4me1", "H3K4me2", "H3K4me3"),
           col=c("#009E73", "black", "orange", "gray"), pch=15, bty="n")

legend("bottomright", inset=c(-0.45,0),
           legend=c("MEP-EryA", "GMP-Granu","GMP-Mono"),
           col="black", lty=1:3, bty="n", lwd=2)


}


dev.off()
warnings()


q()

