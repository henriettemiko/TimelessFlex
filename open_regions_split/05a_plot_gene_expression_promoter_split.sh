#!/bin/bash

##########
#name:          05a_plot_gene_expression_promoter_split.sh
#description:   plots gene expression for promoter clusters split
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 13, 2019
##########


source ../set_variables_hg19.sh

#choose cluster number for promoter clusters here
NUM_CLUSTER=19


SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR_SPLIT/Signal_Generator_promoters
TIMELESS_DIR_PROM=$OPEN_REGIONS_DIR_SPLIT/timeless_promoters

NUM_MARKS=4
NUM_TIME_POINTS=6

MODEL_DIR_PROM=$TIMELESS_DIR_PROM/2_30
NUM_CLUSTER_DIR_PROM=$TIMELESS_DIR_PROM/2_30/$NUM_CLUSTER

mkdir -p $NUM_CLUSTER_DIR_PROM
cd $NUM_CLUSTER_DIR_PROM

TIMES=( D0 D2 D5 D7 D10 D8 )
for TIME in "${TIMES[@]}"
do
    echo $TIME

    RSEM_DIR=$OUTPUT_DIR/RNA/$TIME/RSEM

    RSEM_FILES=($RSEM_DIR/*trimmed.genes.results)
    echo "${RSEM_FILES[@]}"

    for f in "${RSEM_FILES[@]}"
    do
        echo $f

        NAME=$(basename $f .genes.results)

        #remove first line with column headers
        #1st column: gene id, 7th column: expected FPKM
        tail -n +2 $f | cut -f1 > gene_ids.txt
        tail -n +2 $f | cut -f7 > ${NAME}_${TIME}_FPKM.txt

    done

    FPKM_FILES=(*_${TIME}_FPKM.txt)
    echo "${FPKM_FILES[@]}"

    #order of gene ids is the same
    paste "${FPKM_FILES[@]}" > ${TIME}_allFPKM.txt

done

JOIN_FILES=(D*_allFPKM.txt)
echo "${JOIN_FILES[@]}"
#order is D0 D10 D2 D5 D7

#genes ids same for each time point 
paste gene_ids.txt "${JOIN_FILES[@]}" > all_join.txt
#col1 gene id, then for each time point for each replicate FPKM

#until here same for every cluster number

#6 time points, col 36 is cluster assignment
cut -f1,2,3,4,5,6,8,10,36 allCountsNorm_${NUM_CLUSTER}classes.txt > \
    regions_${NUM_CLUSTER}classes.bed 

#for the plot only look at genes that are in exactly one cluster
#it could be that one gene belongs to two or more clusters because it has two 
#or more promoter regions
#enhancer regions not considered because we look at gene expression

#a) unidirectional promoter regions 
# strand is from TSS

awk 'OFS="\t" {if ($7!="." && $8==".") print $0}' \
    regions_${NUM_CLUSTER}classes.bed > regions_${NUM_CLUSTER}classes_uni.txt


#b) bidirectional promoter regions
#two genes assigned, left from + strand, right from - strand
#NEW: for bidirectional take both right and left because we have two TSSs
#write region to col 7 for late

awk 'OFS="\t" {if ($7!="." && $8!=".") 
print $1,$2,$3,$4,$5,"+",$7,".",$9"\n"$1,$2,$3,$4,$5,"-",$8,".",$9}' \
    regions_${NUM_CLUSTER}classes.bed > regions_${NUM_CLUSTER}classes_bidir.txt

cat regions_${NUM_CLUSTER}classes_uni.txt \
    regions_${NUM_CLUSTER}classes_bidir.txt | sort -k1,1 -k2,2n > \
    regions_${NUM_CLUSTER}classes_allprom.txt

#regions are distinct and non-overlapping
#but there could still be regions coming from same gene that are assigned to 
#different clusters

#cut to gene id and cluster number
cut -f 7,9 regions_${NUM_CLUSTER}classes_allprom.txt > \
    regions_${NUM_CLUSTER}classes_allprom2.txt

#genes could be in there multiple times
#multiple lines are only taken once
#gene and class assignment must be same here to be fitting the seen
awk '!seen[$0]++' regions_${NUM_CLUSTER}classes_allprom2.txt > \
    regions_${NUM_CLUSTER}classes_allprom2_nodup.txt

#check if gene is in there and assigned to different clusters
awk 'BEGIN { FS="\t" } { c[$1]++; l[$1,c[$1]]=$0 } END { for (i in c) { 
if (c[i] == 1) for (j = 1; j <= c[i]; j++) print l[i,j] } }' \
    regions_${NUM_CLUSTER}classes_allprom2_nodup.txt | sort > \
    regions_${NUM_CLUSTER}classes_allprom2_nodup_sameassignment.txt
#this is all genes that are unambigously assigned to one clusters, 
#col1 gene id, col2 assignment


join -t $'\t' all_join.txt \
    regions_${NUM_CLUSTER}classes_allprom2_nodup_sameassignment.txt > \
    join_${NUM_CLUSTER}assignments.txt

NUM_TIME_POINTS=5
NUM_REPLICATES=3
#note: number of replicates must be the same for each time point
#for D8 only one replicate, add separately in plotting script

Rscript $SCRIPT_DIR/open_regions_split/plot_gene_expression_cluster_split.r \
    join_${NUM_CLUSTER}assignments.txt ${NUM_TIME_POINTS} ${NUM_REPLICATES} \
    $NUM_CLUSTER $MODEL_DIR_PROM $NUM_CLUSTER_DIR_PROM

#this is all genes that are unambigously assigned to one clusters, 
#col1 gene id, col2 assignment
cp regions_${NUM_CLUSTER}classes_allprom2_nodup_sameassignment.txt \
    genes_assignment.txt


exit

