#!/bin/bash

##########
#name:          05a_plot_gene_expression_promoter_pairs_combined_multi.sh
#description:   plots gene expression for promoters pairs combined multi
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 31, 2019
##########


source ../../set_variables_hg19.sh


NUM_CLUSTERS_PROM_ENH=10
NUM_CLUSTERS_PROM_PROM=5
NUM_CLUSTERS_ENH_ENH=17
NUM_CLUSTER_COMBINED_MULTI=32

TIMELESS_DIR_COMBINED_MULTI=${OPEN_REGIONS_DIR_PAIRS}/timeless_combined_multi

CUR_DIR_COMBINED_MULTI=$TIMELESS_DIR_COMBINED_MULTI/\
${NUM_CLUSTERS_PROM_ENH}_${NUM_CLUSTERS_PROM_PROM}_${NUM_CLUSTERS_ENH_ENH}


cd $CUR_DIR_COMBINED_MULTI


paste allCountsNorm.txt classes-${NUM_CLUSTER_COMBINED_MULTI}_afterEM.txt > \
    allCountsNorm_${NUM_CLUSTER_COMBINED_MULTI}classes.txt


TIMES=( D0 D2 D5 D10 )
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
#order is D0 D10 D2 D5

#genes ids same for each time point 
paste gene_ids.txt "${JOIN_FILES[@]}" > all_join.txt
#col1 gene id, then for each time point for each replicate FPKM

#until here same for every cluster number



#regions are col2-12
cut -f2-7,9,11,56 allCountsNorm_${NUM_CLUSTER_COMBINED_MULTI}classes.txt > \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes1.bed

cut -f29-34,36,38,56 allCountsNorm_${NUM_CLUSTER_COMBINED_MULTI}classes.txt \
    > regions_${NUM_CLUSTER_COMBINED_MULTI}classes2.bed 

#for the plot only look at genes that are in exactly one cluster
#it could be that one gene belongs to two or more clusters because it has two 
#or more promoter regions
#enhancer regions not considered because we look at gene expression

#a) unidirectional promoter regions 
# strand is from TSS

awk 'OFS="\t" {if ($7!="." && $8==".") print $0}' \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes1.bed > \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_uni1.txt

awk 'OFS="\t" {if ($7!="." && $8==".") print $0}' \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes2.bed > \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_uni2.txt

#b) bidirectional promoter regions
#two genes assigned, left from + strand, right from - strand
#NEW: for bidirectional take both right and left because we have two TSSs
#write region to col 7 for late

awk 'OFS="\t" {if ($7!="." && $8!=".") 
print $1,$2,$3,$4,$5,"+",$7,".",$9"\n"$1,$2,$3,$4,$5,"-",$8,".",$9}' \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes1.bed > \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_bidir1.txt

awk 'OFS="\t" {if ($7!="." && $8!=".") 
print $1,$2,$3,$4,$5,"+",$7,".",$9"\n"$1,$2,$3,$4,$5,"-",$8,".",$9}' \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes2.bed > \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_bidir2.txt

cat regions_${NUM_CLUSTER_COMBINED_MULTI}classes_uni1.txt \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_bidir1.txt | \
    sort -k1,1 -k2,2n > \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom1.txt

cat regions_${NUM_CLUSTER_COMBINED_MULTI}classes_uni2.txt \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_bidir2.txt | \
    sort -k1,1 -k2,2n > \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2.txt

#regions are distinct and non-overlapping
#but there could still be regions coming from same gene that are assigned to 
#different clusters

#cut to gene id and cluster number
cut -f 7,9 regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom1.txt > \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_1.txt
cut -f 7,9 regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2.txt > \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_2.txt

#genes could be in there multiple times
#multiple lines are only taken once
#gene and class assignment must be same here to be fitting the seen
awk '!seen[$0]++' regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_1.txt \
    > regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_nodup1.txt
awk '!seen[$0]++' regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_2.txt \
    > regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_nodup2.txt


#check if gene is in there and assigned to different clusters
awk 'BEGIN { FS="\t" } { c[$1]++; l[$1,c[$1]]=$0 } END { for (i in c) { 
if (c[i] == 1) for (j = 1; j <= c[i]; j++) print l[i,j] } }' \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_nodup1.txt | \
    sort > \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_nodup_sameassign1.txt

awk 'BEGIN { FS="\t" } { c[$1]++; l[$1,c[$1]]=$0 } END { for (i in c) { 
if (c[i] == 1) for (j = 1; j <= c[i]; j++) print l[i,j] } }' \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_nodup2.txt | \
    sort > \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_nodup_sameassign2.txt
#this is all genes that are unambigously assigned to one clusters, 
#col1 gene id, col2 assignment


join -t $'\t' all_join.txt \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_nodup_sameassign1.txt \
    > join_${NUM_CLUSTER_COMBINED_MULTI}assignments1.txt

join -t $'\t' all_join.txt \
    regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_nodup_sameassign2.txt \
    > join_${NUM_CLUSTER_COMBINED_MULTI}assignments2.txt

NUM_TIME_POINTS=4
NUM_REPLICATES=3
#note: number of replicates must be the same for each time point

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
plot_gene_expression_cluster_prom-prom.r \
    join_${NUM_CLUSTER_COMBINED_MULTI}assignments1.txt \
    join_${NUM_CLUSTER_COMBINED_MULTI}assignments2.txt ${NUM_TIME_POINTS} \
    ${NUM_REPLICATES} ${NUM_CLUSTER_COMBINED_MULTI} \
    $CUR_DIR_COMBINED_MULTI

#this is all genes that are unambigously assigned to one clusters, 
#col1 gene id, col2 assignment
cp regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_nodup_sameassign1.txt \
    genes_assignment1.txt
cp regions_${NUM_CLUSTER_COMBINED_MULTI}classes_allprom2_nodup_sameassign2.txt \
    genes_assignment2.txt


exit

