
##########
#name:          plot_clusters_coherence_both.r
#description:   plots cluster coherence for counts and fold changes fullset
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          November 21, 2019
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


coherence.counts = read.table(paste0(numcluster.dir, "/", numClusters, "/", type, "_", numClusters, "_coherence_counts.txt"))
coherence.fold = read.table(paste0(numcluster.dir, "/", numClusters, "/", type, "_", numClusters, "_coherence_fold.txt"))


correlations.counts = read.table(paste0(numcluster.dir, "/", numClusters, "/", type, "_", numClusters, "_correlations_counts.txt"))
correlations.fold = read.table(paste0(numcluster.dir, "/", numClusters, "/", type, "_", numClusters, "_correlations_fold.txt"))


print(coherence.counts)
print(coherence.fold)

pdf(paste0(type, "_clusters_", numClusters, "_coherence_both.pdf"), height=6, width=8)
plot(coherence.counts$V1, xlab="clusters", ylab="cluster coherence", ylim=c(0,1), pch=19)
points(coherence.fold$V1, pch=19, col="red")
axis(side=1, at=1:numClusters, labels=1:numClusters)
legend("topleft", legend=c("counts","fold changes"), pch=19, col=c("black", "red"))

dev.off()


print(correlations.counts)
print(correlations.fold)

pdf(paste0(type, "_clusters_", numClusters, "_correlations_both.pdf"), height=6, width=8)
plot(correlations.counts$V1, xlab="clusters", ylab="mean(Person's correlation to cluster mean)", ylim=c(0,1), pch=19)
points(correlations.fold$V1, pch=19, col="red")
axis(side=1, at=1:numClusters, labels=1:numClusters)
legend("topleft", legend=c("counts","fold changes"), pch=19, col=c("black", "red"))

dev.off()
warnings()


q()

