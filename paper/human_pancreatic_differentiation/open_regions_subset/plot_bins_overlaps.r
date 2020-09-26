
##########
#name:          plot_bins_overlaps.r
#description:   plot how many Hi-C bins are overlapping how many peaks
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 12, 2019
##########


library(ggplot2)

setwd(getwd())

args = commandArgs(trailingOnly=TRUE)

t <- read.table(args[1])

print(head(t))

print(min(t$V1))
print(max(t$V1))


pdf("barplot_bins_overlaps.pdf", width=8, height=6)

ggplot(t, aes(x = V1)) + geom_bar() + scale_y_log10() + 
scale_x_continuous(breaks = 1:6, minor_breaks = 1:5) + 
labs(x="number of overlapped open chromatin regions", y = "log10(counts)") + 
ggtitle("How many open chromatin regions are Hi-C bins overlapping?") + 
theme(plot.title = element_text(hjust = 0.5)) + 
theme(legend.title=element_blank())

dev.off()


q()

