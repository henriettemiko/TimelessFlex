#!/bin/bash

##########
#name:          06a_do_ATAC_plot.sh
#description:   plots ATAC signal for clusters
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 11, 2019
##########


source ../set_variables_hg19.sh


#choose cluster number for promoter clusters
NUM_CLUSTER_PROM=18

NUM_MARKS=4
NUM_TIME_POINTS=5

SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR_FULL/Signal_Generator_promoters

TIMELESS_DIR_PROM=$OPEN_REGIONS_DIR_FULL/timeless_promoters
MODEL_DIR_PROM=$TIMELESS_DIR_PROM/2_30
NUM_CLUSTER_DIR_PROM=$TIMELESS_DIR_PROM/2_30/$NUM_CLUSTER_PROM
cd $NUM_CLUSTER_DIR_PROM

ATAC_DIR=$OUTPUT_DIR/ATAC/

#plot ATAC clusters
#for each time point: mean of number of cut sites in region

TIMES=( D0 D2 D5 D7 D10 )
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


NUM_REPLICATES=2
NUM_TIME_POINTS=5 #D0, D2, D5, D7, D10
#we count number of ATAC-seq cut sites in region
#col 11- are cutsites
#col 9 is cluster assignment
paste regions_${NUM_CLUSTER_PROM}_D0_cutsites.bed \
    <( cut -f 11 regions_${NUM_CLUSTER_PROM}_D2_cutsites.bed ) \
    <( cut -f 11  regions_${NUM_CLUSTER_PROM}_D5_cutsites.bed ) \
    <( cut -f 11  regions_${NUM_CLUSTER_PROM}_D7_cutsites.bed ) \
    <( cut -f 11  regions_${NUM_CLUSTER_PROM}_D10_cutsites.bed ) | \
    cut -f 9,11- > regions_${NUM_CLUSTER_PROM}_cutsites_all.bed


############
#plot 4
#to get number of peaks at each time point for the regions in one cluster
#get col7 with merged peaks, split comma to newline, 
#count how often D0, D2, D5, D7,D10 in cluster

for ((i=1; i<=$NUM_CLUSTER_PROM; i++)); do echo $i; \
    awk -v i="$i" -v numcluster="${NUM_CLUSTER_PROM}" 'OFS="\t" {if ($32==i) 
        print $7 > "merged_peaks_regions_"numcluster"_cluster"i".txt"}' \
        $NUM_CLUSTER_DIR_PROM/allCountsNorm_${NUM_CLUSTER_PROM}classes.txt; done

#for all regions together: count how often peaks from D0, D2 etc occur
for ((i=1; i<=$NUM_CLUSTER_PROM; i++)); do echo $i; \
    tr , '\n' < merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}.txt > \
    merged_peaks_separated_regions_${NUM_CLUSTER_PROM}_cluster${i}.txt ; done

for ((i=1; i<=$NUM_CLUSTER_PROM; i++)); do echo $i; \
    cut -d"_" -f 1 \
    merged_peaks_separated_regions_${NUM_CLUSTER_PROM}_cluster${i}.txt | \
    sort | uniq -c > \
    merged_peaks_numbers_regions_${NUM_CLUSTER_PROM}_cluster${i}.txt  ; done

#how to get normalization values
cut -f 7 $NUM_CLUSTER_DIR_PROM/allCountsNorm_${NUM_CLUSTER_PROM}classes.txt | \
    tr , '\n' | cut -d"_" -f 1 | sort | uniq -c > \
    normalization_values_${NUM_CLUSTER_PROM}.txt


##################
#plot 5

for ((i=1; i<=$NUM_CLUSTER_PROM; i++)); do echo $i; tr , '\t' < \
    merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}.txt > \
    merged_peaks_regions_${NUM_CLUSTER_PROM}_for_regions_cluster${i}.txt ; done
#this has peaks in columns now

#for each cluster for each region count how many D0, D2, D5 etc peaks were 
#merged
for ((i=1; i<=$NUM_CLUSTER_PROM; i++)); do echo $i;
    awk -F"D0" '{print NF-1}' \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}.txt > \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}_D0.txt
    awk -F"D2" '{print NF-1}' \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}.txt > \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}_D2.txt
    awk -F"D5" '{print NF-1}' \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}.txt > \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}_D5.txt
    awk -F"D7" '{print NF-1}' \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}.txt > \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}_D7.txt
    awk -F"D10" '{print NF-1}' \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}.txt > \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}_D10.txt

    paste merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}_D0.txt \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}_D2.txt \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}_D5.txt \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}_D7.txt \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}_D10.txt > \
        merged_peaks_regions_${NUM_CLUSTER_PROM}_cluster${i}_all.txt
    #in cols are 0,1,2 etc representing
    # 1 1 1 0 0 means for this region there was one peak from D0, one from D2, 
    #one from D5 and none from D7 or D10 etc

done

Rscript $SCRIPT_DIR/open_regions_fullset/plot_ATAC_clusters.r \
    regions_${NUM_CLUSTER_PROM}_cutsites_all.bed ${NUM_TIME_POINTS} \
    ${NUM_REPLICATES} $NUM_CLUSTER_PROM 


exit

