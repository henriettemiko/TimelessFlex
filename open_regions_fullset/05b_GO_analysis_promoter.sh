#!/bin/bash

##########
#name:          05b_GO_analysis_promoter.sh
#description:   calls GO analysis on promoter clusters
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 10, 2019
##########


source ../set_variables_hg19.sh


#choose cluster number for promoter clusters here
NUM_CLUSTER=14

NUM_MARKS=4
NUM_TIME_POINTS=5
TIMELESS_DIR_PROM=$OPEN_REGIONS_DIR_FULL/timeless_promoters

MODEL_DIR_PROM=$TIMELESS_DIR_PROM/2_30
NUM_CLUSTER_DIR_PROM=$TIMELESS_DIR_PROM/2_30/$NUM_CLUSTER

mkdir -p $NUM_CLUSTER_DIR_PROM
cd $NUM_CLUSTER_DIR_PROM


#get all transcripts from gencode and store gene id in col1 and gene name in 
#col2
awk 'OFS="\t" {if($3=="transcript"){print $10,$16}}' $ANNOTATION | \
    tr -d '";' > gencode_geneid_genename.txt

awk '!seen[$0]++' gencode_geneid_genename.txt > \
    gencode_geneid_genename_nodup.txt

join -t $'\t' <( sort -k1,1 gencode_geneid_genename_nodup.txt ) \
    <( sort -k1,1 genes_assignment.txt) > genes_names_assignment.txt

#split to cluster
for ((i=1; i<=$NUM_CLUSTER; i++)); do echo $i; awk -v i="$i" \
    'OFS="\t" {if ($3==i) print $0 > "genes_names_assignment_class"i".txt"}' \
    genes_names_assignment.txt; done

#split to cluster
for ((i=1; i<=$NUM_CLUSTER; i++)); do echo $i; awk -v i="$i" \
    'OFS="\t" {if ($3==i) print $2 > "names_class"i".txt"}' \
    genes_names_assignment.txt; done

#all genes from all clusters (as background for GO analysis)
cut -f 2 genes_names_assignment.txt > names_all.txt

Rscript $SCRIPT_DIR/open_regions_fullset/go_analysis.r $NUM_CLUSTER \
    $NUM_CLUSTER_DIR_PROM


exit

