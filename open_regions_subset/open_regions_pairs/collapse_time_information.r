
library(plyr)


args = commandArgs(trailingOnly=TRUE)

setwd(getwd())


#collapsing time point information for pairs that occur multiple times in file (one line for one time point before, after one line with all time points comma separated)

t = read.table(args[1])
numClusters=args[2]

#col1 is HiC pair, not important anymore
#col 2-4 and 5-7 are regions of pair
#col 8 is cluster and col 9 time
t2=ddply(t, .(V2,V3,V4,V5,V6,V7,V8), summarize, V9summ=paste(V9,collapse=","))
dim(t2)
print(head(t2))
write.table(t2, file=paste0("pairs_", numClusters, "classes_HiC_collapsed.txt"), sep="\t", quote = F, row.names = F, col.names=F)

t3=ddply(t, .(V2,V3,V4,V5,V6,V7), summarize,  V8summ=paste(V8,collapse=","), V9summ=paste(V9,collapse=","))
dim(t3)
print(head(t3))
write.table(t3, file=paste0("pairs_", numClusters, "classes_HiC_collapsed2.txt"), sep="\t", quote = F, row.names = F, col.names=F)


q()
