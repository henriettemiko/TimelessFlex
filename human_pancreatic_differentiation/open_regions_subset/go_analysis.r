
##########
#name:          go_analysis.r
#description:   GO analysis with gprofiler for promoter clusters subset
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 7, 2019
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
#Ensembl 75: Feb 2014
####

#listMarts(host="aug2014.archive.ensembl.org")
#listMarts(host="feb2014.archive.ensembl.org")

#based on Ensembl Release 75 data
listMarts(host="grch37.ensembl.org")

mart <- useMart("ensembl", host="grch37.ensembl.org", 
                dataset="hsapiens_gene_ensembl")

#p-value strict 0.001
pdf(paste0(numcluster.dir,"/gprofiler_0001_", numClusters, ".pdf"), 
    height=11, width=11)

for (i in 1:numClusters) {

    #read in genes for each cluster (same as for RSEM plots)
    genes <- read.table(paste0(numcluster.dir,"/names_class",i,".txt"))
    print(dim(genes))
    print(head(genes))

    gprofiler(genes, organism = "hsapiens", src_filter = "GO", 
              max_p_value = 1e-3, hier_filtering = "moderate", 
              correction_method = "fdr",include_graph = T, 
              custom_bg=genes.background,
              png_fn=paste0(numcluster.dir,"/gprofiler_0001_", 
                            numClusters, "_", i, ".png"))


    enrichments = gprofiler(genes, organism = "hsapiens", src_filter = "GO", 
                            max_p_value = 1e-3, hier_filtering = "moderate", 
                            correction_method = "fdr", include_graph = F, 
                            custom_bg=genes.background)

    print(dim(enrichments))

    if (dim(enrichments)[1] >0 ){

        enrichments.ordered = enrichments[order(enrichments$domain, 
                                                -enrichments$p.value),]
        enrichments.ordered$ord = 1:nrow(enrichments.ordered)

        gp <- ggplot(data = enrichments.ordered, 
                     mapping = aes(x = reorder(term.name, ord), 
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

    gprofiler(genes, organism = "hsapiens", src_filter = "GO", 
              max_p_value = 1e-2, hier_filtering = "moderate", 
              correction_method = "fdr", include_graph = T, 
              custom_bg=genes.background,
              png_fn=paste0(numcluster.dir,"/gprofiler_001_", 
                            numClusters, "_", i, ".png"))

    enrichments = gprofiler(genes, organism = "hsapiens", src_filter = "GO", 
                            max_p_value = 1e-2, hier_filtering = "moderate", 
                            correction_method = "fdr", include_graph = F, 
                            custom_bg=genes.background)

    print(dim(enrichments))

    if (dim(enrichments)[1] >0 ){

        enrichments.ordered = enrichments[order(enrichments$domain, 
                                                -enrichments$p.value),]
        enrichments.ordered$ord = 1:nrow(enrichments.ordered)

        gp <- ggplot(data = enrichments.ordered, 
                     mapping = aes(x = reorder(term.name, ord), 
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




#p-value strict 0.001 and no hierarchical filtering
pdf(paste0(numcluster.dir,"/gprofiler_0001_", numClusters, "_none.pdf"), 
    height=11, width=11)

for (i in 1:numClusters) {

    #read in genes for each cluster (same as for RSEM plots)
    genes <- read.table(paste0(numcluster.dir,"/names_class",i,".txt"))
    print(dim(genes))
    print(head(genes))

    enrichments = gprofiler(genes, organism = "hsapiens", src_filter = "GO", 
                            max_p_value = 1e-3, hier_filtering = "none", 
                            correction_method = "fdr", include_graph = T, 
                            custom_bg=genes.background,
                            png_fn=paste0(numcluster.dir,"/gprofiler_0001_", 
                                          numClusters,  "_", i, "_none.png"))


    enrichments = gprofiler(genes, organism = "hsapiens", src_filter = "GO", 
                            max_p_value = 1e-3, hier_filtering = "none", 
                            correction_method = "fdr", include_graph = F, 
                            custom_bg=genes.background)
    print(dim(enrichments))

    if (dim(enrichments)[1] >0 ){

        enrichments.ordered = enrichments[order(enrichments$domain, 
                                                -enrichments$p.value),]
        enrichments.ordered$ord = 1:nrow(enrichments.ordered)

        gp <- ggplot(data = enrichments.ordered, 
                     mapping = aes(x = reorder(term.name, ord), 
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




#filtering for set sizes
pdf(paste0(numcluster.dir,"/gprofiler_001_", numClusters, "_filtered.pdf"), 
    height=11, width=11)

for (i in 1:numClusters) {

    #read in genes for each cluster (same as for RSEM plots)
    genes <- read.table(paste0(numcluster.dir,"/names_class",i,".txt"))
    print(dim(genes))
    print(head(genes))

    enrichments = gprofiler(genes, organism = "hsapiens", src_filter = "GO", 
                            max_p_value = 1e-2, hier_filtering = "none", 
                            correction_method = "fdr", include_graph = T, 
                            custom_bg=genes.background,
                            min_set_size=100, max_set_size=2000,
                            png_fn=paste0(numcluster.dir,"/gprofiler_0001_", 
                                          numClusters,  "_", i, "_filtered.png"))


    enrichments = gprofiler(genes, organism = "hsapiens", src_filter = "GO", 
                            max_p_value = 1e-3, hier_filtering = "none", 
                            correction_method = "fdr", include_graph = F, 
                            custom_bg=genes.background)
    print(dim(enrichments))

    if (dim(enrichments)[1] >0 ){

        enrichments.ordered = enrichments[order(enrichments$domain, 
                                                -enrichments$p.value),]
        enrichments.ordered$ord = 1:nrow(enrichments.ordered)

        gp <- ggplot(data = enrichments.ordered, 
                     mapping = aes(x = reorder(term.name, ord), 
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


q()

