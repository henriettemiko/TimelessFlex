#!/bin/bash

##########
#name:          04i_plot_clusters_pairs_initmulti.sh
#description:   calls plotting of clusters for pairs initmulti
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 26, 2019
##########


source ../../set_variables_hg19.sh


NUM_CLUSTERS_PROM_ENH=10
NUM_CLUSTERS_PROM_PROM=4
NUM_CLUSTERS_ENH_ENH=15

TIMELESS_DIR_COMBINED_INITMULTI=${OPEN_REGIONS_DIR_PAIRS}/\
timeless_combined_initmulti


CUR_DIR_COMBINED_INITMULTI=$TIMELESS_DIR_COMBINED_INITMULTI/\
${NUM_CLUSTERS_PROM_ENH}_${NUM_CLUSTERS_PROM_PROM}_${NUM_CLUSTERS_ENH_ENH}

NUM_CLUSTER_COMBINED_INITMULTI=29


NUM_MARKS=4
NUM_TIME_POINTS=4


SIGNAL_GENERATOR_DIR_INIT_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_init_prom-enh

SIGNAL_GENERATOR_DIR_INIT_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_init_prom-prom

SIGNAL_GENERATOR_DIR_INIT_ENH_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_init_enh-enh


SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_prom-enh

SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_prom-prom

SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_enh-enh


cd $CUR_DIR_COMBINED_INITMULTI


cat $SIGNAL_GENERATOR_DIR_INIT_PROM_ENH/allCountsNorm.txt \
    $SIGNAL_GENERATOR_DIR_INIT_PROM_PROM/allCountsNorm.txt \
    $SIGNAL_GENERATOR_DIR_INIT_ENH_ENH/allCountsNorm.txt \
    $SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH/allCountsNorm.txt \
    $SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM/allCountsNorm.txt \
    $SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH/allCountsNorm.txt > allCountsNorm.txt

paste allCountsNorm.txt classes-${NUM_CLUSTER_COMBINED_INITMULTI}_afterEM.txt \
    > allCountsNorm_${NUM_CLUSTER_COMBINED_INITMULTI}classes.txt

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
plot_clusters_pairs.r \
    $NUM_CLUSTER_COMBINED_INITMULTI $NUM_MARKS $NUM_TIME_POINTS \
    $CUR_DIR_COMBINED_INITMULTI $CUR_DIR_COMBINED_INITMULTI \
    $CUR_DIR_COMBINED_INITMULTI "combined_initmulti"


exit

