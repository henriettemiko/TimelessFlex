#!/bin/bash

##########
#name:          05b_GO_analysis_pairs_combined_multi.sh
#description:   calls GO analysis on promoter clusters from pairs combined multi
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 1, 2019
##########


source ../../set_variables_hg19.sh


NUM_CLUSTERS_PROM_ENH=10
NUM_CLUSTERS_PROM_PROM=4
NUM_CLUSTERS_ENH_ENH=15
NUM_CLUSTER_COMBINED_MULTI=29

TIMELESS_DIR_COMBINED_MULTI=${OPEN_REGIONS_DIR_PAIRS}/timeless_combined_multi

CUR_DIR_COMBINED_MULTI=$TIMELESS_DIR_COMBINED_MULTI/\
${NUM_CLUSTERS_PROM_ENH}_${NUM_CLUSTERS_PROM_PROM}_${NUM_CLUSTERS_ENH_ENH}

cd $CUR_DIR_COMBINED_MULTI


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
for ((i=1; i<=$NUM_CLUSTER_COMBINED_MULTI; i++)); do echo $i; awk -v i="$i" \
    'OFS="\t" {if ($3==i) print $0 > "genes_names_assignment_class"i".txt"}' \
    genes_names_assignment.txt; done

#split to cluster
for ((i=1; i<=$NUM_CLUSTER_COMBINED_MULTI; i++)); do echo $i; awk -v i="$i" \
    'OFS="\t" {if ($3==i) print $2 > "names_class"i".txt"}' \
    genes_names_assignment.txt; done

#all genes from all clusters (as background for GO analysis)
cut -f 2 genes_names_assignment.txt > names_all.txt

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/go_analysis.r \
    $NUM_CLUSTER_COMBINED_MULTI \
    $CUR_DIR_COMBINED_MULTI


exit

