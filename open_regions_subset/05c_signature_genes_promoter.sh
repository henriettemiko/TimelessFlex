#!/bin/bash

##########
#name:          05b_signature_genes_promoter.sh
#description:   plots number of signature genes and promoter clusters subset
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 12, 2019
##########


source ../set_variables_hg19.sh


#choose cluster number for promoters here
NUM_CLUSTER=14

NUM_MARKS=4
NUM_TIME_POINTS=4

SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR_SUB/Signal_Generator_promoters
TIMELESS_DIR_PROM=$OPEN_REGIONS_DIR_SUB/timeless_promoters

MODEL_DIR_PROM=$TIMELESS_DIR_PROM/2_30
NUM_CLUSTER_DIR_PROM=$TIMELESS_DIR_PROM/2_30/$NUM_CLUSTER
mkdir -p $NUM_CLUSTER_DIR_PROM
cd $NUM_CLUSTER_DIR_PROM

SIGNATURE_GENES_DIR=$INPUT_DIR/signature_genes
#mkdir -p $SIGNATURE_GENES_DIR
#note: copy signature genes in directory


#############
#take list of signature genes from Xie et al. suppl and see in which clusters 
#they fall
#copy columns from excel sheet into txt file manually

#  685 DE_signature_genes.txt
#  155 GT_signature_genes.txt
#  566 FG_signature_genes.txt (D7 here not used)
#  236 PE_signature_genes.txt
##########


#all genes that are unambiguously(!) assigned and belonging to DE signature 
#list
join -t $'\t' <( cut -f2,3 genes_names_assignment.txt | sort -k1,1) \
    <(sort -k1,1 $SIGNATURE_GENES_DIR/DE_signature_genes.txt) > DE_genes.txt

#to which cluster where these genes assigned, col1 number of genes, 
#col2 assignments
#sum of col1 is DE genes that are unambigously assigned to one cluster
cut -f2 DE_genes.txt | sort -n | uniq -c > DE_genes_clusters.txt

join -t $'\t' <( cut -f2,3 genes_names_assignment.txt | sort -k1,1) \
    <(sort -k1,1 $SIGNATURE_GENES_DIR/GT_signature_genes.txt) > GT_genes.txt

cut -f2 GT_genes.txt | sort -n | uniq -c > GT_genes_clusters.txt

join -t $'\t' <( cut -f2,3 genes_names_assignment.txt | sort -k1,1) \
    <(sort -k1,1 $SIGNATURE_GENES_DIR/PE_signature_genes.txt) > PE_genes.txt

cut -f2 PE_genes.txt | sort -n | uniq -c > PE_genes_clusters.txt

Rscript $SCRIPT_DIR/open_regions_subset/plot_cluster_signature_genes.r \
    $NUM_CLUSTER $NUM_MARKS $NUM_TIME_POINTS $MODEL_DIR_PROM \
    $TIMELESS_DIR_PROM $SIGNAL_GENERATOR_DIR_PROM $NUM_CLUSTER_DIR_PROM


exit

