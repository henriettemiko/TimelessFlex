#!/bin/bash

##########
#name:          06a_plot_ATAC_promoter.sh
#description:   plots ATAC signal for clusters promoters
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          April 19, 2020
##########


source ../set_variables_mm10.sh


#choose cluster number for promoters clusters
NUM_CLUSTER_PROM=16

NUM_MARKS=4
NUM_TIME_POINTS=6

SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR/Signal_Generator_promoters

TIMELESS_DIR_PROM=$OPEN_REGIONS_DIR/timeless_promoters
MODEL_DIR_PROM=$TIMELESS_DIR_PROM/2_30
NUM_CLUSTER_DIR_PROM=$TIMELESS_DIR_PROM/2_30/$NUM_CLUSTER_PROM
cd $NUM_CLUSTER_DIR_PROM

#col36 stores cluster assignments
cut -f1,2,3,4,5,6,8,10,36 allCountsNorm_${NUM_CLUSTER_PROM}classes.txt > \
    regions_${NUM_CLUSTER_PROM}classes.bed 


ATAC_DIR=$OUTPUT_DIR/ATAC/

#plot ATAC clusters
#for each time point: mean of number of cut sites in region

TIMES=( CMP MEP EryA GMP Granu Mono )
for TIME in "${TIMES[@]}"
do
    echo $TIME

    #scaled bedgraphs of 1bp ATAC
    BEDGRAPH_FILES=($ATAC_DIR/${TIME}/bedgraphs/*.bedGraph)
    echo "${BEDGRAPH_FILES[@]}"


    #note: more memory needed: 30G
    #all regions with cluster assignments
    #stores sum of all overlaps of time points, later divide by number of 
    #replicates (here 2) to get mean
    bedtools intersect -c \
        -a $NUM_CLUSTER_DIR_PROM/regions_${NUM_CLUSTER_PROM}classes.bed \
        -b ${BEDGRAPH_FILES[@]} > \
        regions_${NUM_CLUSTER_PROM}_${TIME}_cutsites_unnormalized.bed

done


#col 10 is overlaps
#divide by region length
#col 11 is normalized overlaps
for t in "${TIMES[@]}"; do echo $t; \
    awk -v t="$t" -v numcluster="${NUM_CLUSTER_PROM}" 'OFS="\t" {len=$3-$2; 
        print $0,$10/len > "regions_"numcluster"_"t"_cutsites.bed"}' \
        regions_${NUM_CLUSTER_PROM}_${t}_cutsites_unnormalized.bed; done


NUM_REPLICATES=2 #not for this data set, changes in R script
NUM_TIME_POINTS=4 
#we count number of ATAC-seq cut sites in region
#col 11- are cutsites
#col 9 is cluster assignment
paste regions_${NUM_CLUSTER_PROM}_CMP_cutsites.bed \
    <( cut -f 11 regions_${NUM_CLUSTER_PROM}_MEP_cutsites.bed ) \
    <( cut -f 11  regions_${NUM_CLUSTER_PROM}_EryA_cutsites.bed ) \
    <( cut -f 11  regions_${NUM_CLUSTER_PROM}_GMP_cutsites.bed ) \
    <( cut -f 11  regions_${NUM_CLUSTER_PROM}_Granu_cutsites.bed ) \
    <( cut -f 11  regions_${NUM_CLUSTER_PROM}_Mono_cutsites.bed ) | \
    cut -f 9,11- > regions_${NUM_CLUSTER_PROM}_cutsites_all.bed


Rscript $SCRIPT_DIR/open_regions/plot_ATAC_clusters.r \
    regions_${NUM_CLUSTER_PROM}_cutsites_all.bed ${NUM_TIME_POINTS} \
    ${NUM_REPLICATES} $NUM_CLUSTER_PROM 


exit

