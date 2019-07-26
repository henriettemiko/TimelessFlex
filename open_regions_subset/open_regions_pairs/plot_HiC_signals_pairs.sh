#!/bin/bash

##########
#name:          plot_HiC_signals_pairs.sh
#description:   plots HiC signal for clusters for pairs
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 26, 2019
##########


NUM_CLUSTER=$1
SCRIPT_DIR=$2
OUTPUT_DIR=$3
SIGNAL_GENERATOR_DIR=$4
MODEL_DIR=$5
OPEN_REGIONS_DIR=$6
TYPE=$7

echo $TYPE

FILENAME=""
if [[ $TYPE =~ "multi" ]]
then
    FILENAME="_afterEM"
fi

echo $FILENAME


paste $SIGNAL_GENERATOR_DIR/allCountsNorm.txt \
    $MODEL_DIR/\
classes-${NUM_CLUSTER}${FILENAME}.txt > \
    allCountsNorm_${NUM_CLUSTER}classes.txt

#here cut regions from pairs
#keep peak information from ATAC col 8 and cluster assignment col 56
#cut -f2-7,8,9,11,56 allCountsNorm_${NUM_CLUSTER}classes.txt > \
#    prom_regions_${NUM_CLUSTER}classes.bed

#col1 is Hi-C pair information
#col2-4 is region1 information
#col 29-31 is region2
#col56 is cluster assignment

cut -f1,2-4,29-31,56 allCountsNorm_${NUM_CLUSTER}classes.txt > \
    pairs_${NUM_CLUSTER}classes.bed


paste <(cut -d"." -f2- ${OPEN_REGIONS_DIR}/mypairs_all_notunique.txt) <(cut -d"." -f1 ${OPEN_REGIONS_DIR}/mypairs_all_notunique.txt) > mypairs_all_notunique_time.txt


#get Hi-C time information




join -t $'\t' -1 1 -2 1 <(sort -k1,1 pairs_${NUM_CLUSTER}classes.bed ) <(sort -k1,1 mypairs_all_notunique_time.txt ) > pairs_${NUM_CLUSTER}classes_HiC.txt




#collape time information with R
Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/collapse_time_information.r pairs_${NUM_CLUSTER}classes_HiC.txt $NUM_CLUSTER
#col 7 is cluster assignment for pair and col 8 time points from Hi-C interactions




#split pairs in separate files for each cluster
for ((i=1; i<=$NUM_CLUSTER; i++)); do echo $i; awk -v i="$i" -v numcluster="${NUM_CLUSTER}" 'OFS="\t" {if ($7==i) print $8 > "HiC_"numcluster"_cluster"i".txt"}' pairs_${NUM_CLUSTER}classes_HiC_collapsed.txt; done



#replace , with newline
for ((i=1; i<=$NUM_CLUSTER; i++)); do echo $i; tr , '\n' < HiC_${NUM_CLUSTER}_cluster${i}.txt > HiC_${NUM_CLUSTER}_cluster${i}_newline.txt ; done

#counts how often D0, D2, D5, D10 occur in each cluster
for ((i=1; i<=$NUM_CLUSTER; i++)); do echo $i; sort HiC_${NUM_CLUSTER}_cluster${i}_newline.txt | uniq -c > HiC_${NUM_CLUSTER}_cluster${i}_newline_counts.txt ; done


#how to get normalization values
#for all unique prom-enh pairs, how many D0,D2,D5,D10?
cut -f 8 pairs_${NUM_CLUSTER}classes_HiC_collapsed.txt | tr , '\n' | sort | uniq -c > normalization_values_${NUM_CLUSTER}.txt


#if I want values and not one number

for ((i=1; i<=$NUM_CLUSTER; i++)); do echo $i;


awk -F"D0" '{print NF-1}' HiC_${NUM_CLUSTER}_cluster${i}.txt > HiC_${NUM_CLUSTER}_cluster${i}_D0.txt
awk -F"D2" '{print NF-1}' HiC_${NUM_CLUSTER}_cluster${i}.txt > HiC_${NUM_CLUSTER}_cluster${i}_D2.txt
awk -F"D5" '{print NF-1}' HiC_${NUM_CLUSTER}_cluster${i}.txt > HiC_${NUM_CLUSTER}_cluster${i}_D5.txt
awk -F"D10" '{print NF-1}' HiC_${NUM_CLUSTER}_cluster${i}.txt > HiC_${NUM_CLUSTER}_cluster${i}_D10.txt


paste HiC_${NUM_CLUSTER}_cluster${i}_D0.txt HiC_${NUM_CLUSTER}_cluster${i}_D2.txt HiC_${NUM_CLUSTER}_cluster${i}_D5.txt HiC_${NUM_CLUSTER}_cluster${i}_D10.txt > HiC_${NUM_CLUSTER}_cluster${i}_all.txt 


#matrix
# 0 1 0 1 #interactions in D2 and D10, not in D0 and D5

done


#plot in R
Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_HiC_clusters.r HiC_${NUM_CLUSTER}_cluster${i}_newline_counts.txt HiC_${NUM_CLUSTER}_cluster${i}_all.txt normalization_values_${NUM_CLUSTER}.txt ${NUM_CLUSTER}


exit


