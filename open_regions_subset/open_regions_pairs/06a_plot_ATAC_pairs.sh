#!/bin/bash

##########
#name:          06a_plot_ATAC_pairs.sh
#description:   plots ATAC signal for clusters for pairs
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 24, 2019
##########


source ../../set_variables_hg19.sh


#choose cluster number for multi prom-enh pairs here
NUM_CLUSTER_MULTI_PROM_ENH=10

#choose cluster number for multi prom-prom pairs here
NUM_CLUSTER_MULTI_PROM_PROM=5

#choose cluster number for multi enh-enh pairs here
NUM_CLUSTER_MULTI_ENH_ENH=12


SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}_noordering/\
Signal_Generator_multi_prom-enh
TIMELESS_DIR_MULTI_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}_noordering/timeless_multi_prom-enh

MODEL_DIR_MULTI_PROM_ENH=$TIMELESS_DIR_MULTI_PROM_ENH/\
$NUM_CLUSTER_MULTI_PROM_ENH
NUM_CLUSTER_DIR_MULTI_PROM_ENH=$TIMELESS_DIR_MULTI_PROM_ENH/\
$NUM_CLUSTER_MULTI_PROM_ENH
mkdir -p $NUM_CLUSTER_DIR_MULTI_PROM_ENH

cd $NUM_CLUSTER_DIR_MULTI_PROM_ENH

paste $SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH/allCountsNorm.txt \
    $MODEL_DIR_MULTI_PROM_ENH/\
classes-${NUM_CLUSTER_MULTI_PROM_ENH}_afterEM.txt > \
    allCountsNorm_${NUM_CLUSTER_MULTI_PROM_ENH}classes.txt

#here cut regions from pairs
#keep peak information from ATAC col 8 and cluster assignment col 56
cut -f2-7,8,9,11,56 allCountsNorm_${NUM_CLUSTER_MULTI_PROM_ENH}classes.txt > \
    prom_regions_${NUM_CLUSTER_MULTI_PROM_ENH}classes.bed

cut -f29-34,35,36,38,56 allCountsNorm_${NUM_CLUSTER_MULTI_PROM_ENH}classes.txt > \
    enh_regions_${NUM_CLUSTER_MULTI_PROM_ENH}classes.bed 


ATAC_DIR=$OUTPUT_DIR/ATAC/

#plot ATAC clusters
#for each time point: mean of number of cut sites in region

TIMES=( D0 D2 D5 D10 )
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

#    bedtools intersect -c \
#        -a prom_regions_${NUM_CLUSTER_MULTI_PROM_ENH}classes.bed \
#        -b ${BEDGRAPH_FILES[@]} > \
#        prom_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_${TIME}_cutsites_unnormalized.bed

#    bedtools intersect -c \
#        -a enh_regions_${NUM_CLUSTER_MULTI_PROM_ENH}classes.bed \
#        -b ${BEDGRAPH_FILES[@]} > \
#        enh_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_${TIME}_cutsites_unnormalized.bed

done


#col 11 is overlaps
#divide by region length
#col 12 is normalized overlaps
for t in "${TIMES[@]}"; do echo $t; \
    awk -v t="$t" -v numcluster="${NUM_CLUSTER_MULTI_PROM_ENH}" 'OFS="\t" {len=$3-$2; 
print $0,$11/len > "prom_regions_"numcluster"_"t"_cutsites.bed"}' \
    prom_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_${t}_cutsites_unnormalized.bed; done

for t in "${TIMES[@]}"; do echo $t; \
    awk -v t="$t" -v numcluster="${NUM_CLUSTER_MULTI_PROM_ENH}" 'OFS="\t" {len=$3-$2; 
print $0,$11/len > "enh_regions_"numcluster"_"t"_cutsites.bed"}' \
    enh_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_${t}_cutsites_unnormalized.bed; done


NUM_REPLICATES=2
NUM_TIME_POINTS=4 #D0, D2, D5, D10
#we count number of ATAC-seq cut sites in region
#col 12- are cutsites
#col 10 is cluster assignment
paste prom_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_D0_cutsites.bed \
    <( cut -f 12 prom_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_D2_cutsites.bed ) \
    <( cut -f 12  prom_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_D5_cutsites.bed ) \
    <( cut -f 12  prom_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_D10_cutsites.bed ) | \
    cut -f 10,12- > prom_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cutsites_all.bed


paste enh_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_D0_cutsites.bed \
    <( cut -f 12 enh_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_D2_cutsites.bed ) \
    <( cut -f 12  enh_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_D5_cutsites.bed ) \
    <( cut -f 12  enh_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_D10_cutsites.bed ) | \
    cut -f 10,12- > enh_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cutsites_all.bed

############
#plot 4
#to get number of peaks at each time point for the regions in one cluster
#get col7 with merged peaks, split comma to newline, 
#count how often D0, D2, D5, D10 in cluster

for ((i=1; i<=$NUM_CLUSTER_MULTI_PROM_ENH; i++)); do
    awk -v i="$i" -v numcluster="${NUM_CLUSTER_MULTI_PROM_ENH}" 'OFS="\t" {if ($10==i) 
print $7 > "prom_merged_peaks_regions_"numcluster"_cluster"i".txt"}' \
    prom_regions_${NUM_CLUSTER_MULTI_PROM_ENH}classes.bed; done
#allCountsNorm_${NUM_CLUSTER_MULTI_PROM_ENH}classes.txt; done



for ((i=1; i<=$NUM_CLUSTER_MULTI_PROM_ENH; i++)); do
    awk -v i="$i" -v numcluster="${NUM_CLUSTER_MULTI_PROM_ENH}" 'OFS="\t" {if ($10==i) 
print $7 > "enh_merged_peaks_regions_"numcluster"_cluster"i".txt"}' \
    enh_regions_${NUM_CLUSTER_MULTI_PROM_ENH}classes.bed; done


#for all regions together: count how often peaks from D0, D2 etc occur
for ((i=1; i<=$NUM_CLUSTER_MULTI_PROM_ENH; i++)); do
    tr , '\n' < prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt > \
    prom_merged_peaks_separated_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt ; done

for ((i=1; i<=$NUM_CLUSTER_MULTI_PROM_ENH; i++)); do
    tr , '\n' < enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt > \
    enh_merged_peaks_separated_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt ; done

for ((i=1; i<=$NUM_CLUSTER_MULTI_PROM_ENH; i++)); do
    cut -d"_" -f 1 \
    prom_merged_peaks_separated_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt | \
    sort | uniq -c > \
    prom_merged_peaks_numbers_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt  ; done

for ((i=1; i<=$NUM_CLUSTER_MULTI_PROM_ENH; i++)); do
    cut -d"_" -f 1 \
    enh_merged_peaks_separated_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt | \
    sort | uniq -c > \
    enh_merged_peaks_numbers_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt  ; done

#how to get normalization values
cut -f 7 prom_regions_${NUM_CLUSTER_MULTI_PROM_ENH}classes.bed | \
    tr , '\n' | cut -d"_" -f 1 | sort | uniq -c > \
    prom_normalization_values_${NUM_CLUSTER_MULTI_PROM_ENH}.txt

cut -f 7 enh_regions_${NUM_CLUSTER_MULTI_PROM_ENH}classes.bed | \
    tr , '\n' | cut -d"_" -f 1 | sort | uniq -c > \
    enh_normalization_values_${NUM_CLUSTER_MULTI_PROM_ENH}.txt


##################
#plot 5


#promoter
for ((i=1; i<=$NUM_CLUSTER_MULTI_PROM_ENH; i++)); do tr , '\t' < \
    prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt > \
    prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_for_regions_cluster${i}.txt ; done
#this has peaks in columns now

#for each cluster for each region count how many D0, D2, D5 etc peaks were 
#merged
for ((i=1; i<=$NUM_CLUSTER_MULTI_PROM_ENH; i++)); do
    awk -F"D0" '{print NF-1}' \
        prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt > \
        prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D0.txt
    awk -F"D2" '{print NF-1}' \
        prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt > \
        prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D2.txt
    awk -F"D5" '{print NF-1}' \
        prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt > \
        prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D5.txt
    awk -F"D10" '{print NF-1}' \
        prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt > \
        prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D10.txt

    paste prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D0.txt \
        prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D2.txt \
        prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D5.txt \
        prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D10.txt > \
        prom_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_all.txt
    #in cols are 0,1,2 etc representing
    # 1 1 1 0 means for this region there was one peak from D0, one from D2, 
    #one from D5 and none from D10 etc
done

#enhancer
for ((i=1; i<=$NUM_CLUSTER_MULTI_PROM_ENH; i++)); do tr , '\t' < \
    enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt > \
    enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_for_regions_cluster${i}.txt ; done
#this has peaks in columns now

#for each cluster for each region count how many D0, D2, D5 etc peaks were 
#merged
for ((i=1; i<=$NUM_CLUSTER_MULTI_PROM_ENH; i++)); do
    awk -F"D0" '{print NF-1}' \
        enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt > \
        enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D0.txt
    awk -F"D2" '{print NF-1}' \
        enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt > \
        enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D2.txt
    awk -F"D5" '{print NF-1}' \
        enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt > \
        enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D5.txt
    awk -F"D10" '{print NF-1}' \
        enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}.txt > \
        enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D10.txt

    paste enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D0.txt \
        enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D2.txt \
        enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D5.txt \
        enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_D10.txt > \
        enh_merged_peaks_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cluster${i}_all.txt
    #in cols are 0,1,2 etc representing
    # 1 1 1 0 means for this region there was one peak from D0, one from D2, 
    #one from D5 and none from D10 etc
done

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_ATAC_clusters_pairs.r \
    prom_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cutsites_all.bed \
    enh_regions_${NUM_CLUSTER_MULTI_PROM_ENH}_cutsites_all.bed ${NUM_TIME_POINTS} \
    ${NUM_REPLICATES} $NUM_CLUSTER_MULTI_PROM_ENH "multi_prom-enh"


exit

