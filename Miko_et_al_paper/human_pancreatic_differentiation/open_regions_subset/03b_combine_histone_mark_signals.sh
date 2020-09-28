#!/bin/bash

##########
#name:          03b_combine_histone_mark_signals.sh
#description:   compute histone mark signals at feature regions subset
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 12, 2019
##########


source ../set_variables_hg19.sh


FEATURE_REGIONS_DIR=$OUTPUT_DIR/open_regions_subset/get_feature_regions
PROM_REGIONS=$FEATURE_REGIONS_DIR/all_ext_500_final_promoters.bed
ENH_REGIONS=$FEATURE_REGIONS_DIR/all_ext_500_final_enhancers.bed

echo "prom regions files: $PROM_REGIONS"
echo `wc -l $PROM_REGIONS`
echo "enh regions files: $ENH_REGIONS"
echo `wc -l $ENH_REGIONS`


SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR_SUB/Signal_Generator_promoters
SIGNAL_GENERATOR_DIR_ENH=$OPEN_REGIONS_DIR_SUB/Signal_Generator_enhancers

NUM_MARKS=4
NUM_TIME_POINTS=5

cd $SIGNAL_GENERATOR_DIR_PROM

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

MAX_FILES_PROM=($SIGNAL_GENERATOR_DIR_PROM/*/*/max_counts.txt)
echo "${MAX_FILES_PROM[@]}"

paste <(sort-bed-typical $PROM_REGIONS) "${MAX_FILES_PROM[@]}" > \
    all_max_counts_prom.txt

Rscript $SCRIPT_DIR/open_regions_subset/get_normalized_FC_nofilter.r \
    all_max_counts_prom.txt "_prom" "promoters"

#the files allCountsNorm and allFold are ordered in D0, D2, D5, D10 and 
#D2/D0, D5/D2, D10/D5
#stored in SIGNAL_GENERATOR_DIR
#order of histone marks as before: H3K27ac H3K27me3 H3K4me1 H3K4me3

cd $SIGNAL_GENERATOR_DIR_ENH

MAX_FILES_ENH=($SIGNAL_GENERATOR_DIR_ENH/*/*/max_counts.txt)
echo "${MAX_FILES_ENH[@]}"

paste <(sort-bed-typical $ENH_REGIONS) "${MAX_FILES_ENH[@]}" > \
    all_max_counts_enh.txt

Rscript $SCRIPT_DIR/open_regions_subset/get_normalized_FC_nofilter.r \
    all_max_counts_enh.txt "_enh" "enhancers"


exit

