#!/bin/bash

##########
#name:          05b_signature_genes_promoter_multi.sh
#description:   plots number of signature genes and promoter clusters pairs multi
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 30, 2019
##########


source ../../set_variables_hg19.sh


#choose cluster number for multi prom-enh pairs here
NUM_CLUSTER_MULTI_PROM_ENH=10

#choose cluster number for multi prom-prom pairs here
NUM_CLUSTER_MULTI_PROM_PROM=5


NUM_MARKS=4
NUM_TIME_POINTS=4



#multi prom-enh

SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_prom-enh
TIMELESS_DIR_MULTI_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
timeless_multi_prom-enh

MODEL_DIR_MULTI_PROM_ENH=$TIMELESS_DIR_MULTI_PROM_ENH/\
$NUM_CLUSTER_MULTI_PROM_ENH
NUM_CLUSTER_DIR_MULTI_PROM_ENH=$TIMELESS_DIR_MULTI_PROM_ENH/\
$NUM_CLUSTER_MULTI_PROM_ENH
mkdir -p $NUM_CLUSTER_DIR_MULTI_PROM_ENH

cd $NUM_CLUSTER_DIR_MULTI_PROM_ENH


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

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
plot_cluster_signature_genes_pairs.r \
    $NUM_CLUSTER_MULTI_PROM_ENH $NUM_MARKS $NUM_TIME_POINTS \
    $MODEL_DIR_MULTI_PROM_ENH \
    $SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH \
    $NUM_CLUSTER_DIR_MULTI_PROM_ENH "multi_prom-enh"



#multi prom-prom

SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_prom-prom
TIMELESS_DIR_MULTI_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/\
timeless_multi_prom-prom

MODEL_DIR_MULTI_PROM_PROM=$TIMELESS_DIR_MULTI_PROM_PROM/\
$NUM_CLUSTER_MULTI_PROM_PROM
NUM_CLUSTER_DIR_MULTI_PROM_PROM=$TIMELESS_DIR_MULTI_PROM_PROM/\
$NUM_CLUSTER_MULTI_PROM_PROM
mkdir -p $NUM_CLUSTER_DIR_MULTI_PROM_PROM

cd $NUM_CLUSTER_DIR_MULTI_PROM_PROM


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


#keep sides separately

join -t $'\t' <( sort -k1,1 gencode_geneid_genename_nodup.txt ) \
    <( sort -k1,1 genes_assignment1.txt) > genes_names_assignment1.txt

join -t $'\t' <( sort -k1,1 gencode_geneid_genename_nodup.txt ) \
    <( sort -k1,1 genes_assignment2.txt) > genes_names_assignment2.txt


#together

join -t $'\t' <( sort -k1,1 gencode_geneid_genename_nodup.txt ) \
    <( sort -k1,1 genes_assignment.txt) > genes_names_assignment.txt


#all genes that are unambiguously(!) assigned and belonging to DE signature 
#list
join -t $'\t' <( cut -f2,3 genes_names_assignment.txt | sort -k1,1) \
    <(sort -k1,1 $SIGNATURE_GENES_DIR/DE_signature_genes.txt) > DE_genes.txt

join -t $'\t' <( cut -f2,3 genes_names_assignment1.txt | sort -k1,1) \
    <(sort -k1,1 $SIGNATURE_GENES_DIR/DE_signature_genes.txt) > DE_genes1.txt

join -t $'\t' <( cut -f2,3 genes_names_assignment2.txt | sort -k1,1) \
    <(sort -k1,1 $SIGNATURE_GENES_DIR/DE_signature_genes.txt) > DE_genes2.txt

#to which cluster where these genes assigned, col1 number of genes, 
#col2 assignments
#sum of col1 is DE genes that are unambigously assigned to one cluster
cut -f2 DE_genes1.txt | sort -n | uniq -c > DE_genes_clusters1.txt
cut -f2 DE_genes2.txt | sort -n | uniq -c > DE_genes_clusters2.txt


join -t $'\t' <( cut -f2,3 genes_names_assignment.txt | sort -k1,1) \
    <(sort -k1,1 $SIGNATURE_GENES_DIR/GT_signature_genes.txt) > GT_genes.txt
join -t $'\t' <( cut -f2,3 genes_names_assignment1.txt | sort -k1,1) \
    <(sort -k1,1 $SIGNATURE_GENES_DIR/GT_signature_genes.txt) > GT_genes1.txt
join -t $'\t' <( cut -f2,3 genes_names_assignment2.txt | sort -k1,1) \
    <(sort -k1,1 $SIGNATURE_GENES_DIR/GT_signature_genes.txt) > GT_genes2.txt

cut -f2 GT_genes1.txt | sort -n | uniq -c > GT_genes_clusters1.txt
cut -f2 GT_genes2.txt | sort -n | uniq -c > GT_genes_clusters2.txt


join -t $'\t' <( cut -f2,3 genes_names_assignment.txt | sort -k1,1) \
    <(sort -k1,1 $SIGNATURE_GENES_DIR/PE_signature_genes.txt) > PE_genes.txt
join -t $'\t' <( cut -f2,3 genes_names_assignment1.txt | sort -k1,1) \
    <(sort -k1,1 $SIGNATURE_GENES_DIR/PE_signature_genes.txt) > PE_genes1.txt
join -t $'\t' <( cut -f2,3 genes_names_assignment2.txt | sort -k1,1) \
    <(sort -k1,1 $SIGNATURE_GENES_DIR/PE_signature_genes.txt) > PE_genes2.txt

cut -f2 PE_genes1.txt | sort -n | uniq -c > PE_genes_clusters1.txt
cut -f2 PE_genes2.txt | sort -n | uniq -c > PE_genes_clusters2.txt

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
plot_cluster_signature_genes_pairs_prom-prom.r \
    $NUM_CLUSTER_MULTI_PROM_PROM $NUM_MARKS $NUM_TIME_POINTS \
    $MODEL_DIR_MULTI_PROM_PROM \
    $SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM \
    $NUM_CLUSTER_DIR_MULTI_PROM_PROM "multi_prom-prom"


exit

