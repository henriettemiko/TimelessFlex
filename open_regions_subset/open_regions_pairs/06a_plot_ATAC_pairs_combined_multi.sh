#!/bin/bash

##########
#name:          06a_plot_ATAC_pairs_combined_multi.sh
#description:   plots ATAC signal for clusters for pairs combined multi
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 31, 2019
##########


source ../../set_variables_hg19.sh



NUM_CLUSTERS_PROM_ENH=10
NUM_CLUSTERS_PROM_PROM=4
NUM_CLUSTERS_ENH_ENH=15

TIMELESS_DIR_COMBINED_MULTI=${OPEN_REGIONS_DIR_PAIRS}/timeless_combined_multi


CUR_DIR_COMBINED_MULTI=$TIMELESS_DIR_COMBINED_MULTI/\
${NUM_CLUSTERS_PROM_ENH}_${NUM_CLUSTERS_PROM_PROM}_${NUM_CLUSTERS_ENH_ENH}

NUM_CLUSTER_COMBINED_MULTI=32

cd $CUR_DIR_COMBINED_MULTI

#combined multi

$SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_ATAC_signals_pairs.sh \
    $NUM_CLUSTER_COMBINED_MULTI $SCRIPT_DIR $OUTPUT_DIR \
    $CUR_DIR_COMBINED_MULTI $CUR_DIR_COMBINED_MULTI \
    "combined_multi"
    

exit

