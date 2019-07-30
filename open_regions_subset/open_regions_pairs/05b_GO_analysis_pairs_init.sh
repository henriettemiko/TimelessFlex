#!/bin/bash

##########
#name:          05b_GO_analysis_pairs_init.sh
#description:   calls GO analysis on promoter clusters from pairs init
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 24, 2019
##########


source ../../set_variables_hg19.sh



#choose cluster number for init prom-enh pairs here
NUM_CLUSTER_INIT_PROM_ENH=17

#choose cluster number for init prom-prom pairs here
NUM_CLUSTER_INIT_PROM_PROM=7



#init prom-enh

SIGNAL_GENERATOR_DIR_INIT_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_init_prom-enh
TIMELESS_DIR_INIT_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
timeless_init_prom-enh

MODEL_DIR_INIT_PROM_ENH=$TIMELESS_DIR_INIT_PROM_ENH/2_30/
NUM_CLUSTER_DIR_INIT_PROM_ENH=$TIMELESS_DIR_INIT_PROM_ENH/\
$NUM_CLUSTER_INIT_PROM_ENH
mkdir -p $NUM_CLUSTER_DIR_INIT_PROM_ENH

cd $NUM_CLUSTER_DIR_INIT_PROM_ENH


#get all transcripts from gencode and store gene id in col1 and gene name in 
#col2
awk 'OFS="\t" {if($3=="transcript"){print $10,$16}}' $ANNOTATION | \
    tr -d '";' > gencode_geneid_genename.txt

awk '!seen[$0]++' gencode_geneid_genename.txt > \
    gencode_geneid_genename_nodup.txt

join -t $'\t' <( sort -k1,1 gencode_geneid_genename_nodup.txt ) \
    <( sort -k1,1 genes_assignment.txt) > genes_names_assignment.txt

#split to cluster
for ((i=1; i<=$NUM_CLUSTER_INIT_PROM_ENH; i++)); do echo $i; awk -v i="$i" \
    'OFS="\t" {if ($3==i) print $0 > "genes_names_assignment_class"i".txt"}' \
    genes_names_assignment.txt; done

#split to cluster
for ((i=1; i<=$NUM_CLUSTER_INIT_PROM_ENH; i++)); do echo $i; awk -v i="$i" \
    'OFS="\t" {if ($3==i) print $2 > "names_class"i".txt"}' \
    genes_names_assignment.txt; done

#all genes from all clusters (as background for GO analysis)
cut -f 2 genes_names_assignment.txt > names_all.txt

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/go_analysis.r \
    $NUM_CLUSTER_INIT_PROM_ENH \
    $NUM_CLUSTER_DIR_INIT_PROM_ENH



#init prom-prom


SIGNAL_GENERATOR_DIR_INIT_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_init_prom-prom
TIMELESS_DIR_INIT_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/\
timeless_init_prom-prom

MODEL_DIR_INIT_PROM_PROM=$TIMELESS_DIR_INIT_PROM_PROM/2_30/
NUM_CLUSTER_DIR_INIT_PROM_PROM=$TIMELESS_DIR_INIT_PROM_PROM/\
$NUM_CLUSTER_INIT_PROM_PROM
mkdir -p $NUM_CLUSTER_DIR_INIT_PROM_PROM

cd $NUM_CLUSTER_DIR_INIT_PROM_PROM


#get all transcripts from gencode and store gene id in col1 and gene name in 
#col2
awk 'OFS="\t" {if($3=="transcript"){print $10,$16}}' $ANNOTATION | \
    tr -d '";' > gencode_geneid_genename.txt

awk '!seen[$0]++' gencode_geneid_genename.txt > \
    gencode_geneid_genename_nodup.txt


#concatenate gene assignments from both sides
cat genes_assignment1.txt genes_assignment2.txt > genes_assignment.txt



join -t $'\t' <( sort -k1,1 gencode_geneid_genename_nodup.txt ) \
    <( sort -k1,1 genes_assignment.txt) > genes_names_assignment.txt

#split to cluster
for ((i=1; i<=$NUM_CLUSTER_INIT_PROM_PROM; i++)); do echo $i; awk -v i="$i" \
    'OFS="\t" {if ($3==i) print $0 > "genes_names_assignment_class"i".txt"}' \
    genes_names_assignment.txt; done

#split to cluster
for ((i=1; i<=$NUM_CLUSTER_INIT_PROM_PROM; i++)); do echo $i; awk -v i="$i" \
    'OFS="\t" {if ($3==i) print $2 > "names_class"i".txt"}' \
    genes_names_assignment.txt; done

#all genes from all clusters (as background for GO analysis)
cut -f 2 genes_names_assignment.txt > names_all.txt

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/go_analysis.r \
    $NUM_CLUSTER_INIT_PROM_PROM \
    $NUM_CLUSTER_DIR_INIT_PROM_PROM


exit

