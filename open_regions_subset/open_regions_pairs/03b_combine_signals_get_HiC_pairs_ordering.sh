#!/bin/bash

##########
#name:          03b_combine_signals_get_HiC_pairs_ordering.sh
#description:   get Hi-C pairs and order prom-prom and enh-enh pairs
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 23, 2019
##########


source ../../set_variables_hg19.sh


FEATURE_REGIONS_DIR=$OUTPUT_DIR/open_regions_subset/get_feature_regions
PROM_REGIONS=$FEATURE_REGIONS_DIR/all_ext_500_final_promoters.bed
ENH_REGIONS=$FEATURE_REGIONS_DIR/all_ext_500_final_enhancers.bed

echo "prom regions files: $PROM_REGIONS"
echo `wc -l $PROM_REGIONS`
echo "enh regions files: $ENH_REGIONS"
echo `wc -l $ENH_REGIONS`

SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR_SUB/Signal_Generator_promoters
SIGNAL_GENERATOR_DIR_ENH=$OPEN_REGIONS_DIR_SUB/Signal_Generator_enhancers

#OPEN_REGIONS_DIR_PAIRS=$OPEN_REGIONS_DIR_SUB/open_regions_pairs
#mkdir -p $OPEN_REGIONS_DIR_PAIRS

cd $OPEN_REGIONS_DIR_PAIRS

MAX_FILES_PROM=($SIGNAL_GENERATOR_DIR_PROM/*/*/max_counts.txt)

paste <(sort-bed-typical $PROM_REGIONS) "${MAX_FILES_PROM[@]}" > \
    all_max_counts_prom.txt


MAX_FILES_ENH=($SIGNAL_GENERATOR_DIR_ENH/*/*/max_counts.txt)

paste <(sort-bed-typical $ENH_REGIONS) "${MAX_FILES_ENH[@]}" > \
    all_max_counts_enh.txt


#feature regions with counts
cat all_max_counts_prom.txt all_max_counts_enh.txt | sort -k1,1 -k2,2n > \
    all_max_counts_prom_enh.txt


#combine counts from Signal Generator for each mark and for each time point 
#in defined order:
#H3K27ac D0 D10 D2 D5
#H3K27me3 D0 D10 D2 D5
#H3K4me1 D0 D10 D2 D5
#H3K4me3 D0 D10 D2 D5

#ordering of time points: D0 D2 D5 D10
#ordering of columns: col1 col3 col4 col2
#FC of time points D2/D0 D5/D2 D10/D5
#FC columns: col3/col1 col4/col3 col2/col4


###information about pairs of Hi-C bins that are interacting###

#3 sets of Hi-C bin pairs:
#pairs where both bins overlap at least one open chromatin region (ALL pairs)
#pairs where both bins overlap exactely one open chromatin region 
#(1:1 relationship, INIT pairs)
#pairs where one or both bins overlap more than one open chromatin regions 
#(MULTI pairs)


#ALL PAIRS:
#take all bins with overlap (one or multiple), but take each bin once, 
#then get pair infomation
cut -f 4 "$OPEN_REGIONS_DIR_SUB/combine_peaks/"\
"all_bins_unique_overlapping_peaks.bed" | sort | uniq | tr , "\n" | \
    sed s/.right//g | sed s/.left//g | sort | uniq -d > \
    mypairs_all_notunique.txt
#ALL PAIRS FROM UNCOLLAPSED BINS WHERE BIN OVERLAPS AT LEAST ONE OPEN 
#CHROMATIN REGION FULLY
#here we still have time information, so not unique

#cut time information
cut -d"." -f2- mypairs_all_notunique.txt | sort | uniq > mypairs_all_unique.txt

#INIT PAIRS:
#in file bins can be there multiple times because they can have multiple 
#overlaps with open chromatin regions
#take unique bins that have one overlap with open chromatin region (uniq -u)
#col 4 is bin with pair information(chr*.*.*.left/right*)
#bin need to be there twice (once .left and once .right)
#if it is only there once the interaction bin is lost

cut -f 4 "$OPEN_REGIONS_DIR_SUB/combine_peaks/"\
"all_bins_unique_overlapping_peaks.bed" | sort | uniq -u | tr , "\n" | \
    sed s/.right//g | sed s/.left//g | sort | uniq -d > \
    mypairs_init_notunique.txt
#INIT SET

cut -d"." -f2- mypairs_init_notunique.txt | sort | uniq > \
    mypairs_init_unique.txt


#MULTI PAIRS:
#ALL PAIRS minus INIT PAIRS
comm -23 <( sort mypairs_all_notunique.txt ) \
    <( sort mypairs_init_notunique.txt ) > mypairs_multi_notunique.txt
#MULTI SET


comm -23 <( sort mypairs_all_unique.txt ) <( sort mypairs_init_unique.txt ) > \
    mypairs_multi_unique.txt


#############

#we need bins to get from Hi-C to feature region information

#here we have peak identifier, not feature region (but is same)
#identifier needed for join
awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$7"."$8,$10,$11,$12}' \
    "$OPEN_REGIONS_DIR_SUB/combine_peaks/\
"all_bins_notunique_overlapping_peaks.bed > \
    all_bins_notunique_overlapping_peaks_identifier.bed



#combine region counts to pairs again before clustering pairs
####INIT PAIRS

awk 'OFS="\t" {print $0".left\n"$0".right"}' mypairs_init_notunique.txt > \
    mypairs_init_leftright.txt
#here we have Hi-C informations

#col 4 is comma separated list for unique bins
#so better take not unique bins here
#here join on column time.chromosome.start.end.left/right
join -t $'\t' -1 1 -2 4 \
    -o 2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 \
    <(sort mypairs_init_leftright.txt) \
    <(sort -k4,4 all_bins_notunique_overlapping_peaks_identifier.bed) > \
    mypairs_init_leftright_bins.bed
#here only 13 cols because peaks and not feature regions!


join -t $'\t' -1 4 -2 10 \
    -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,\
1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,2.4 \
<( sed s/.ext//g all_max_counts_prom_enh.txt | sort -k4,4 ) \
<(grep left mypairs_init_leftright_bins.bed | sed s/.left//g | sort -k10,10 ) \
> max_counts_init_left.txt

join -t $'\t' -1 4 -2 10 \
    -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,\
1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,2.4 \
    <( sed s/.ext//g all_max_counts_prom_enh.txt | sort -k4,4 ) \
    <(grep right mypairs_init_leftright_bins.bed | sed s/.right//g | \
    sort -k10,10 ) > max_counts_init_right.txt


join -t $'\t' -1 28 -2 28  <(sort -k28,28 max_counts_init_left.txt ) \
    <( sort -k28,28 max_counts_init_right.txt) > max_counts_init_pairs.txt
#not unique


#get unique pairs by cutting part of col1 (time information)
cut -d "." -f 2- max_counts_init_pairs.txt | sort| uniq > \
    max_counts_init_pairs_unique.txt

#HiC info is col 1, region left is 2-12, counts left 13-28, 
#region right 29-39, counts right 40-55
cut -f 13-28,39-54 max_counts_init_pairs_unique.txt > \
    max_counts_init_pairs_unique_counts.txt


###INIT###
#split into prom-enh, prom-prom and enh-enh

#prom-prom
awk 'OFS="\t" {if ($9!="." && $36!=".") print $0}' \
    max_counts_init_pairs_unique.txt > \
    max_counts_init_pairs_unique_prom-prom.txt

###
#reorder prom-prom pairs: right side has higher K27ac value at D10

awk 'OFS="\t" {if ($16<=$43) {print $0}}' \
    max_counts_init_pairs_unique_prom-prom.txt > \
    max_counts_init_pairs_unique_prom-prom_correctorder_pt1.txt
awk 'OFS="\t" {if ($16>$43) {print $0}}' \
    max_counts_init_pairs_unique_prom-prom.txt > \
    max_counts_init_pairs_unique_prom-prom_wrongorder.txt

cut -f 2-28 max_counts_init_pairs_unique_prom-prom_wrongorder.txt > \
    left_prom-prom.txt
cut -f 1,29- max_counts_init_pairs_unique_prom-prom_wrongorder.txt > \
    right_prom-prom.txt
paste right_prom-prom.txt left_prom-prom.txt > \
    max_counts_init_pairs_unique_prom-prom_correctorder_pt2.txt

cat max_counts_init_pairs_unique_prom-prom_correctorder_pt1.txt \
    max_counts_init_pairs_unique_prom-prom_correctorder_pt2.txt | \
    sort -k1,1 -k2,2n > \
    max_counts_init_pairs_unique_prom-prom.txt

###

cut -f 13-28,40-55 max_counts_init_pairs_unique_prom-prom.txt > \
    max_counts_init_pairs_unique_prom-prom_counts.txt


#enh-enh
awk 'OFS="\t" {if ($9=="."  && $36=="." ) print $0}' \
    max_counts_init_pairs_unique.txt > max_counts_init_pairs_unique_enh-enh.txt

###
#reorder enh-enh pairs: right side has higher K27ac value at D10

awk 'OFS="\t" {if ($16<=$43) {print $0}}' \
    max_counts_init_pairs_unique_enh-enh.txt > \
    max_counts_init_pairs_unique_enh-enh_correctorder_pt1.txt
awk 'OFS="\t" {if ($16>$43) {print $0}}' \
    max_counts_init_pairs_unique_enh-enh.txt > \
    max_counts_init_pairs_unique_enh-enh_wrongorder.txt

cut -f 2-28 max_counts_init_pairs_unique_enh-enh_wrongorder.txt > \
    left_enh-enh.txt
cut -f 1,29- max_counts_init_pairs_unique_enh-enh_wrongorder.txt > \
    right_enh-enh.txt
paste right_enh-enh.txt left_enh-enh.txt > \
    max_counts_init_pairs_unique_enh-enh_correctorder_pt2.txt

cat max_counts_init_pairs_unique_enh-enh_correctorder_pt1.txt \
    max_counts_init_pairs_unique_enh-enh_correctorder_pt2.txt | \
    sort -k1,1 -k2,2n > \
    max_counts_init_pairs_unique_enh-enh.txt

###

cut -f 13-28,40-55 max_counts_init_pairs_unique_enh-enh.txt > \
    max_counts_init_pairs_unique_enh-enh_counts.txt


#prom-enh
awk 'OFS="\t" {if ($9=="." && $36!="." ) print $0}' \
    max_counts_init_pairs_unique.txt > \
    max_counts_init_pairs_unique_enh-prom_pt1.txt

awk 'OFS="\t" {if ($9!="." && $36=="." ) print $0}' \
    max_counts_init_pairs_unique.txt > \
    max_counts_init_pairs_unique_prom-enh_pt1.txt

#change order from enh-prom to prom-enh
cut -f 2-28 max_counts_init_pairs_unique_enh-prom_pt1.txt > left_init.txt
cut -f 1,29- max_counts_init_pairs_unique_enh-prom_pt1.txt > right_init.txt
paste right_init.txt left_init.txt > \
    max_counts_init_pairs_unique_prom-enh_pt2.txt

cat max_counts_init_pairs_unique_prom-enh_pt1.txt \
    max_counts_init_pairs_unique_prom-enh_pt2.txt | sort -k1,1 -k2,2n > \
    max_counts_init_pairs_unique_prom-enh.txt 

cut -f 13-28,40-55 max_counts_init_pairs_unique_prom-enh.txt > \
    max_counts_init_pairs_unique_prom-enh_counts.txt


#compute normalized counts and fold changes

SIGNAL_GENERATOR_DIR_INIT_PROM_ENH=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_init_prom-enh
mkdir -p $SIGNAL_GENERATOR_DIR_INIT_PROM_ENH
cd $SIGNAL_GENERATOR_DIR_INIT_PROM_ENH

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
get_normalized_FC_nofilter_pairs.r \
$OPEN_REGIONS_DIR_PAIRS/max_counts_init_pairs_unique_prom-enh.txt \
"init_prom-enh"


SIGNAL_GENERATOR_DIR_INIT_PROM_PROM=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_init_prom-prom
mkdir -p $SIGNAL_GENERATOR_DIR_INIT_PROM_PROM
cd $SIGNAL_GENERATOR_DIR_INIT_PROM_PROM

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
get_normalized_FC_nofilter_pairs.r \
    $OPEN_REGIONS_DIR_PAIRS/max_counts_init_pairs_unique_prom-prom.txt \
    "init_prom-prom"


SIGNAL_GENERATOR_DIR_INIT_ENH_ENH=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_init_enh-enh
mkdir -p $SIGNAL_GENERATOR_DIR_INIT_ENH_ENH
cd $SIGNAL_GENERATOR_DIR_INIT_ENH_ENH

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
get_normalized_FC_nofilter_pairs.r \
    $OPEN_REGIONS_DIR_PAIRS/max_counts_init_pairs_unique_enh-enh.txt \
    "init_enh-enh"


###MULTI###

cd $OPEN_REGIONS_DIR_PAIRS

awk 'OFS="\t" {print $0".left\n"$0".right"}' mypairs_multi_notunique.txt > \
    mypairs_multi_leftright.txt

join -t $'\t' -1 1 -2 4 \
    -o 2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 \
    <(sort mypairs_multi_leftright.txt) \
    <(sort -k4,4 all_bins_notunique_overlapping_peaks_identifier.bed) > \
    mypairs_multi_leftright_bins.bed


#here join on peak identifier chromosome.start.ext
join -t $'\t' -1 4 -2 10 \
    -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,\
1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,2.4  \
<( sed s/.ext//g all_max_counts_prom_enh.txt | sort -k4,4 ) \
<(grep left mypairs_multi_leftright_bins.bed | sed s/.left//g | \
sort -k10,10 ) > max_counts_multi_left.txt

join -t $'\t' -1 4 -2 10 \
    -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,\
1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,2.4 \
<( sed s/.ext//g all_max_counts_prom_enh.txt | sort -k4,4 ) \
<(grep right mypairs_multi_leftright_bins.bed | sed s/.right//g | \
sort -k10,10 ) > max_counts_multi_right.txt

join -t $'\t' -1 28 -2 28  <(sort -k28,28 max_counts_multi_left.txt ) \
    <( sort -k28,28 max_counts_multi_right.txt) > max_counts_multi_pairs.txt 
#not unique


#get unique pairs by cutting part of col1 (time information)
cut -d "." -f 2- max_counts_multi_pairs.txt | sort | uniq > \
    max_counts_multi_pairs_unique.txt

cut -f 13-28,39-54 max_counts_multi_pairs_unique.txt > \
    max_counts_multi_pairs_unique_counts.txt


#split into multi prom-enh, prom-prom and enh-enh

#prom-prom
awk 'OFS="\t" {if ($9!="." && $36!=".") print $0}' \
    max_counts_multi_pairs_unique.txt > \
    max_counts_multi_pairs_unique_prom-prom.txt

cut -f 13-28,40-55 max_counts_multi_pairs_unique_prom-prom.txt > \
    max_counts_multi_pairs_unique_prom-prom_counts.txt


#enh-enh
awk 'OFS="\t" {if ($9=="."  && $36=="." ) print $0}' \
    max_counts_multi_pairs_unique.txt > \
    max_counts_multi_pairs_unique_enh-enh.txt

cut -f 13-28,40-55 max_counts_multi_pairs_unique_enh-enh.txt > \
    max_counts_multi_pairs_unique_enh-enh_counts.txt


#prom-enh
awk 'OFS="\t" {if ($9=="." && $36!="." ) print $0}' \
    max_counts_multi_pairs_unique.txt > \
    max_counts_multi_pairs_unique_enh-prom_pt1.txt

awk 'OFS="\t" {if ($9!="." && $36=="." ) print $0}' \
    max_counts_multi_pairs_unique.txt > \
    max_counts_multi_pairs_unique_prom-enh_pt1.txt

cut -f 2-28 max_counts_multi_pairs_unique_enh-prom_pt1.txt > left_multi.txt
cut -f 1,29- max_counts_multi_pairs_unique_enh-prom_pt1.txt > right_multi.txt
paste right_multi.txt left_multi.txt > \
    max_counts_multi_pairs_unique_prom-enh_pt2.txt


cat max_counts_multi_pairs_unique_prom-enh_pt1.txt \
    max_counts_multi_pairs_unique_prom-enh_pt2.txt | sort -k1,1 -k2,2n > \
    max_counts_multi_pairs_unique_prom-enh.txt 

cut -f 13-28,40-55 max_counts_multi_pairs_unique_prom-enh.txt > \
    max_counts_multi_pairs_unique_prom-enh_counts.txt


#compute normalized counts and fold changes

SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_multi_prom-enh
mkdir -p $SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH
cd $SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
get_normalized_FC_nofilter_pairs.r \
$OPEN_REGIONS_DIR_PAIRS/max_counts_multi_pairs_unique_prom-enh.txt \
"multi_prom-enh"


SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_multi_prom-prom
mkdir -p $SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM
cd $SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
get_normalized_FC_nofilter_pairs.r \
$OPEN_REGIONS_DIR_PAIRS/max_counts_multi_pairs_unique_prom-prom.txt \
"multi_prom-prom"


SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_multi_enh-enh
mkdir -p $SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH
cd $SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
get_normalized_FC_nofilter_pairs.r \
$OPEN_REGIONS_DIR_PAIRS/max_counts_multi_pairs_unique_enh-enh.txt \
"multi_enh-enh"


exit






###########NOT NEEDED ANYMORE############


###information of Hi-C bins belonging to Hi-C pairs###

grep -f mypairs_all_notunique.txt \
    $OPEN_REGIONS_DIR_SUB/combine_peaks/all_bins_unique_overlapping_peaks.bed \
    > all_bins_merged_overlapping_peaks_mypairs_all.bed
#not unique

grep -f mypairs_all_unique.txt \
    $OPEN_REGIONS_DIR_SUB/combine_peaks/all_bins_unique_overlapping_peaks.bed \
    > all_bins_merged_overlapping_peaks_mypairs_all_unique.bed

grep -f mypairs_init.txt \
    $OPEN_REGIONS_DIR_SUB/combine_peaks/all_bins_unique_overlapping_peaks.bed \
    > all_bins_merged_overlapping_peaks_mypairs_init.bed
#unique

grep -f mypairs_multi.txt \
    $OPEN_REGIONS_DIR_SUB/combine_peaks/all_bins_unique_overlapping_peaks.bed \
    > all_bins_merged_overlapping_peaks_mypairs_multi.bed


##get unique bins belonging to pairs

cut -f 1-4 all_bins_merged_overlapping_peaks_mypairs_all.bed | sort | uniq > \
    all_bins_merged_peaks_mypairs_all_unique.bed

cut -f 1-4 all_bins_merged_overlapping_peaks_mypairs_init.bed | sort | uniq > \
    all_bins_merged_peaks_mypairs_init_unique.bed

cut -f 1-4 all_bins_merged_overlapping_peaks_mypairs_multi.bed | sort | uniq \
    > all_bins_merged_peaks_mypairs_multi_unique.bed


###information about open chromatin regions that are overlapped by Hi-C bins###


#col 6: HiC time information
cut -f 10 all_bins_merged_overlapping_peaks_mypairs_init.bed | sort | uniq > \
    mypairs_init_peaks.txt

cut -f 10 all_bins_merged_overlapping_peaks_mypairs_all.bed | sort | uniq > \
    mypairs_all_peaks.txt

cut -f 10 all_bins_merged_overlapping_peaks_mypairs_multi.bed | sort | uniq > \
    mypairs_multi_peaks.txt
#all open chromatin regions that are overlapped by bins from pairs with 
#multiple overlaps
#these are unique peaks belonging to the pairs


exit

