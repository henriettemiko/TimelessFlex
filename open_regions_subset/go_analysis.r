
##########
#name:          go_analysis.r
#description:   GO analysis with gprofiler for promoter clusters subset
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 12, 2019
##########


library("biomaRt")
library(gProfileR)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

numClusters = as.numeric(args[1])
numcluster.dir=args[2]

#as background all genes from all clusters (but only the unambigously assigned
#genes that were also used for RSEM plots)
genes.background <- scan(paste0(numcluster.dir,"/names_all.txt"), 
                         character(), sep="\n")

print(dim(genes.background))

####
#https://biit.cs.ut.ee/gprofiler/page/archives/
#Human hg19 = GRCh37 = Ensembl versions 59-75
#https://www.ensembl.org/info/website/archives/index.html
#Ensembl 76: Aug 2014
####

listMarts(host="aug2014.archive.ensembl.org")

mart <- useMart("ensembl", host="aug2014.archive.ensembl.org", 
                dataset="hsapiens_gene_ensembl")

#p-value strict 0.001
pdf(paste0(numcluster.dir,"/gprofiler_0001_", numClusters, ".pdf"), 
    height=11, width=11)

for (i in 1:numClusters) {

    #read in genes for each cluster (same as for RSEM plots)
    genes <- read.table(paste0(numcluster.dir,"/names_class",i,".txt"))
    print(dim(genes))
    print(head(genes))

    enrichments = gprofiler(genes, organism = "hsapiens", src_filter = "GO", 
                            max_p_value = 1e-3, hier_filtering = "moderate", 
                            correction_method = "fdr", include_graph = F, 
                            custom_bg=genes.background)

    print(dim(enrichments))

    if (dim(enrichments)[1] >0 ){

        enrichments.ordered = enrichments[order(enrichments$domain, 
                                                -enrichments$p.value),]
        enrichments.ordered$orden = 1:nrow(enrichments.ordered)

        gp <- ggplot(data = enrichments.ordered, 
                     mapping = aes(x = reorder(term.name, orden), 
                                   y = -log10(p.value), 
                                   fill = factor(domain))) +
geom_bar(stat = "identity") + coord_flip() + 
scale_x_discrete(drop = F) + xlab(label = "") + 
ylab("-log10(p-value)") + 
ggtitle(paste0("Enrichment analysis for ", nrow(genes), 
               " genes from cluster ", i,
               "\np-value of enrichment < 0.001")) + 
theme(legend.title = element_blank())

        print(gp)

    }

}
dev.off()


#p-value 0.01
pdf(paste0(numcluster.dir,"/gprofiler_001_", numClusters, ".pdf"), 
    height=11, width=11)

for (i in 1:numClusters) {

    genes <- read.table(paste0(numcluster.dir,"/names_class",i,".txt"))

    print(dim(genes))
    print(head(genes))

    enrichments = gprofiler(genes, organism = "hsapiens", src_filter = "GO", 
                            max_p_value = 1e-2, hier_filtering = "moderate", 
                            correction_method = "fdr", include_graph = F, 
                            custom_bg=genes.background)

    print(dim(enrichments))

    if (dim(enrichments)[1] >0 ){

        enrichments.ordered = enrichments[order(enrichments$domain, 
                                                -enrichments$p.value),]
        enrichments.ordered$orden = 1:nrow(enrichments.ordered)

        gp <- ggplot(data = enrichments.ordered, 
                     mapping = aes(x = reorder(term.name, orden), 
                                   y = -log10(p.value), 
                                   fill = factor(domain))) +
geom_bar(stat = "identity") + 
coord_flip() + scale_x_discrete(drop = F) + xlab(label = "") + 
ylab("-log10(p-value)") + 
ggtitle(paste0("Enrichment analysis for ", nrow(genes), 
               " genes from cluster ", i,
               "\np-value of enrichment < 0.01")) + 
theme(legend.title = element_blank())

        print(gp)

    }

}

dev.off()


q()

