#!/bin/bash

##########
#name:          04g_plot_clusters_pairs_multi.sh
#description:   calls plotting of clusters for pairs multi
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 26, 2019
##########


source ../../set_variables_hg19.sh


NUM_CLUSTERS_PROM_ENH=10
NUM_CLUSTERS_PROM_PROM=4
NUM_CLUSTERS_ENH_ENH=15

TIMELESS_DIR_COMBINED_MULTI=${OPEN_REGIONS_DIR_PAIRS}/timeless_combined_multi


CUR_DIR_COMBINED_MULTI=$TIMELESS_DIR_COMBINED_MULTI/${NUM_CLUSTERS_PROM_ENH}_${NUM_CLUSTERS_PROM_PROM}_${NUM_CLUSTERS_ENH_ENH}

NUM_CLUSTER_COMBINED_MULTI=29


NUM_MARKS=4
NUM_TIME_POINTS=4

TIMELESS_DIR_COMBINED_MULTI=${OPEN_REGIONS_DIR_PAIRS}/timeless_combined_multi

SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_prom-enh

SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_prom-prom

SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_enh-enh





cd $CUR_DIR_COMBINED_MULTI


cat $SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH/allCountsNorm.txt $SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM/allCountsNorm.txt $SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH/allCountsNorm.txt > allCountsNorm.txt

paste allCountsNorm.txt classes-${NUM_CLUSTER_COMBINED_MULTI}_afterEM.txt > \
    allCountsNorm_${NUM_CLUSTER_COMBINED_MULTI}classes.txt

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
plot_clusters_pairs.r \
    $NUM_CLUSTER_COMBINED_MULTI $NUM_MARKS $NUM_TIME_POINTS \
    $CUR_DIR_COMBINED_MULTI $CUR_DIR_COMBINED_MULTI \
    $CUR_DIR_COMBINED_MULTI "combined_multi"


exit

