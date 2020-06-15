#!/bin/bash

##########
#name:          06b_plot_HiC_pairs_combined_initmulti.sh
#description:   plots Hi-C signal for clusters for pairs combined initmulti
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 31, 2019
##########


source ../../set_variables_hg19.sh


NUM_CLUSTERS_PROM_ENH=10
NUM_CLUSTERS_PROM_PROM=5
NUM_CLUSTERS_ENH_ENH=17
NUM_CLUSTER_COMBINED_INITMULTI=32

TIMELESS_DIR_COMBINED_INITMULTI=${OPEN_REGIONS_DIR_PAIRS}/timeless_combined_initmulti

CUR_DIR_COMBINED_INITMULTI=$TIMELESS_DIR_COMBINED_INITMULTI/\
${NUM_CLUSTERS_PROM_ENH}_${NUM_CLUSTERS_PROM_PROM}_${NUM_CLUSTERS_ENH_ENH}


cd $CUR_DIR_COMBINED_INITMULTI

$SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_HiC_signals_pairs.sh \
    $NUM_CLUSTER_COMBINED_INITMULTI $SCRIPT_DIR $OUTPUT_DIR \
    $CUR_DIR_COMBINED_INITMULTI $CUR_DIR_COMBINED_INITMULTI \
    ${OPEN_REGIONS_DIR_PAIRS} \
    "combined_initmulti"
    

exit

