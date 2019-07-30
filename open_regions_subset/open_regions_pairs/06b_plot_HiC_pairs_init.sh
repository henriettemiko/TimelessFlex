#!/bin/bash

##########
#name:          06b_plot_HiC_pairs_init.sh
#description:   plots Hi-C signal for clusters for pairs init
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 26, 2019
##########


source ../../set_variables_hg19.sh


#choose cluster number for init prom-enh pairs here
NUM_CLUSTER_INIT_PROM_ENH=17

#choose cluster number for init prom-prom pairs here
NUM_CLUSTER_INIT_PROM_PROM=7

#choose cluster number for init enh-enh pairs here
NUM_CLUSTER_INIT_ENH_ENH=23



#init prom-enh

SIGNAL_GENERATOR_DIR_INIT_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_init_prom-enh
TIMELESS_DIR_INIT_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/timeless_init_prom-enh

MODEL_DIR_INIT_PROM_ENH=$TIMELESS_DIR_INIT_PROM_ENH/2_30/
NUM_CLUSTER_DIR_INIT_PROM_ENH=$TIMELESS_DIR_INIT_PROM_ENH/\
$NUM_CLUSTER_INIT_PROM_ENH
mkdir -p $NUM_CLUSTER_DIR_INIT_PROM_ENH

cd $NUM_CLUSTER_DIR_INIT_PROM_ENH

$SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_HiC_signals_pairs.sh \
    $NUM_CLUSTER_INIT_PROM_ENH $SCRIPT_DIR $OUTPUT_DIR \
    $SIGNAL_GENERATOR_DIR_INIT_PROM_ENH $MODEL_DIR_INIT_PROM_ENH \
    ${OPEN_REGIONS_DIR_PAIRS} \
    "init_prom-enh"
    


#init prom-prom

SIGNAL_GENERATOR_DIR_INIT_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_init_prom-prom
TIMELESS_DIR_INIT_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/timeless_init_prom-prom

MODEL_DIR_INIT_PROM_PROM=$TIMELESS_DIR_INIT_PROM_PROM/2_30/
NUM_CLUSTER_DIR_INIT_PROM_PROM=$TIMELESS_DIR_INIT_PROM_PROM/\
$NUM_CLUSTER_INIT_PROM_PROM
mkdir -p $NUM_CLUSTER_DIR_INIT_PROM_PROM

cd $NUM_CLUSTER_DIR_INIT_PROM_PROM

$SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_HiC_signals_pairs.sh \
    $NUM_CLUSTER_INIT_PROM_PROM $SCRIPT_DIR $OUTPUT_DIR \
    $SIGNAL_GENERATOR_DIR_INIT_PROM_PROM $MODEL_DIR_INIT_PROM_PROM \
    ${OPEN_REGIONS_DIR_PAIRS} \
    "init_prom-prom"
    


#init enh-enh

SIGNAL_GENERATOR_DIR_INIT_ENH_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_init_enh-enh
TIMELESS_DIR_INIT_ENH_ENH=${OPEN_REGIONS_DIR_PAIRS}/timeless_init_enh-enh

MODEL_DIR_INIT_ENH_ENH=$TIMELESS_DIR_INIT_ENH_ENH/2_30/
NUM_CLUSTER_DIR_INIT_ENH_ENH=$TIMELESS_DIR_INIT_ENH_ENH/\
$NUM_CLUSTER_INIT_ENH_ENH
mkdir -p $NUM_CLUSTER_DIR_INIT_ENH_ENH

cd $NUM_CLUSTER_DIR_INIT_ENH_ENH

$SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_HiC_signals_pairs.sh \
    $NUM_CLUSTER_INIT_ENH_ENH $SCRIPT_DIR $OUTPUT_DIR \
    $SIGNAL_GENERATOR_DIR_INIT_ENH_ENH $MODEL_DIR_INIT_ENH_ENH \
    ${OPEN_REGIONS_DIR_PAIRS} \
    "init_enh-enh"
    

exit

