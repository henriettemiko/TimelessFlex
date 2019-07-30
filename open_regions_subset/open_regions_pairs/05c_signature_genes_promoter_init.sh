#!/bin/bash

##########
#name:          05b_signature_genes_promoter.sh
#description:   plots number of signature genes and promoter clusters pairs
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 24, 2019
##########


source ../../set_variables_hg19.sh


#choose cluster number for init prom-enh pairs here
NUM_CLUSTER_INIT_PROM_ENH=10

#choose cluster number for init prom-prom pairs here
NUM_CLUSTER_INIT_PROM_PROM=5


NUM_MARKS=4
NUM_TIME_POINTS=4


SIGNAL_GENERATOR_DIR_INIT_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_init_prom-enh
TIMELESS_DIR_INIT_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
timeless_init_prom-enh

MODEL_DIR_INIT_PROM_ENH=$TIMELESS_DIR_INIT_PROM_ENH/2_30/
NUM_CLUSTER_DIR_INIT_PROM_ENH=$TIMELESS_DIR_INIT_PROM_ENH/2_30/\
$NUM_CLUSTER_INIT_PROM_ENH
mkdir -p $NUM_CLUSTER_DIR_INIT_PROM_ENH

cd $NUM_CLUSTER_DIR_INIT_PROM_ENH


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

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_cluster_signature_genes_pairs.r \
    $NUM_CLUSTER_INIT_PROM_ENH $NUM_MARKS $NUM_TIME_POINTS \
    $TIMELESS_DIR_INIT_PROM_ENH $SIGNAL_GENERATOR_DIR_INIT_PROM_ENH \
    $NUM_CLUSTER_DIR_INIT_PROM_ENH "init_prom-enh"


exit
