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
NUM_CLUSTER_MULTI_PROM_PROM=4

#choose cluster number for multi enh-enh pairs here
NUM_CLUSTER_MULTI_ENH_ENH=15



#multi prom-enh

SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_prom-enh
TIMELESS_DIR_MULTI_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/timeless_multi_prom-enh

MODEL_DIR_MULTI_PROM_ENH=$TIMELESS_DIR_MULTI_PROM_ENH/\
$NUM_CLUSTER_MULTI_PROM_ENH
NUM_CLUSTER_DIR_MULTI_PROM_ENH=$TIMELESS_DIR_MULTI_PROM_ENH/\
$NUM_CLUSTER_MULTI_PROM_ENH
mkdir -p $NUM_CLUSTER_DIR_MULTI_PROM_ENH

cd $NUM_CLUSTER_DIR_MULTI_PROM_ENH

$SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_ATAC_signals_pairs.sh \
    $NUM_CLUSTER_MULTI_PROM_ENH $SCRIPT_DIR $OUTPUT_DIR \
    $SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH $MODEL_DIR_MULTI_PROM_ENH \
    "multi_prom-enh"
    



#multi prom-prom

SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_prom-prom
TIMELESS_DIR_MULTI_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/timeless_multi_prom-prom

MODEL_DIR_MULTI_PROM_PROM=$TIMELESS_DIR_MULTI_PROM_PROM/\
$NUM_CLUSTER_MULTI_PROM_PROM
NUM_CLUSTER_DIR_MULTI_PROM_PROM=$TIMELESS_DIR_MULTI_PROM_PROM/\
$NUM_CLUSTER_MULTI_PROM_PROM
mkdir -p $NUM_CLUSTER_DIR_MULTI_PROM_PROM

cd $NUM_CLUSTER_DIR_MULTI_PROM_PROM

$SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_ATAC_signals_pairs.sh \
    $NUM_CLUSTER_MULTI_PROM_PROM $SCRIPT_DIR $OUTPUT_DIR \
    $SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM $MODEL_DIR_MULTI_PROM_PROM \
    "multi_prom-prom"
    


#multi enh-enh

SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_enh-enh
TIMELESS_DIR_MULTI_ENH_ENH=${OPEN_REGIONS_DIR_PAIRS}/timeless_multi_enh-enh

MODEL_DIR_MULTI_ENH_ENH=$TIMELESS_DIR_MULTI_ENH_ENH/\
$NUM_CLUSTER_MULTI_ENH_ENH
NUM_CLUSTER_DIR_MULTI_ENH_ENH=$TIMELESS_DIR_MULTI_ENH_ENH/\
$NUM_CLUSTER_MULTI_ENH_ENH
mkdir -p $NUM_CLUSTER_DIR_MULTI_ENH_ENH

cd $NUM_CLUSTER_DIR_MULTI_ENH_ENH

$SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_ATAC_signals_pairs.sh \
    $NUM_CLUSTER_MULTI_ENH_ENH $SCRIPT_DIR $OUTPUT_DIR \
    $SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH $MODEL_DIR_MULTI_ENH_ENH \
    "multi_enh-enh"
    



exit
