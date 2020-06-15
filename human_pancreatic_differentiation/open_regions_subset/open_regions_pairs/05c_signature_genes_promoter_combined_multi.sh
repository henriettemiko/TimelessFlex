#!/bin/bash

##########
#name:          05b_signature_genes_promoter_combined_multi.sh
#description:   plots number of signature genes for pairs combined multi
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


NUM_MARKS=4
NUM_TIME_POINTS=4


cd $CUR_DIR_COMBINED_MULTI


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
    $NUM_CLUSTER_COMBINED_MULTI $NUM_MARKS $NUM_TIME_POINTS \
    $CUR_DIR_COMBINED_MULTI \
    $CUR_DIR_COMBINED_MULTI \
    $CUR_DIR_COMBINED_MULTI "combined_multi"


exit

