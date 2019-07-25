#!/bin/bash

##########
#name:          04f_cluster_multi_with_combined.sh
#description:   calls clustering scripts for pairs multi with combined model
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 25, 2019
##########


source ../../set_variables_hg19.sh

#choose numbers here
NUM_CLUSTERS_PROM_ENH=10
NUM_CLUSTERS_PROM_PROM=5
NUM_CLUSTERS_ENH_ENH=12


NUM_MARKS=4
NUM_FC=6

SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_prom-enh

SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_prom-prom

SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH=${OPEN_REGIONS_DIR_PAIRS}/\
Signal_Generator_multi_enh-enh


TIMELESS_DIR_MULTI_PROM_ENH=${OPEN_REGIONS_DIR_PAIRS}/timeless_multi_prom-enh
TIMELESS_DIR_MULTI_PROM_PROM=${OPEN_REGIONS_DIR_PAIRS}/timeless_multi_prom-prom
TIMELESS_DIR_MULTI_ENH_ENH=${OPEN_REGIONS_DIR_PAIRS}/timeless_multi_enh-enh


TIMELESS_DIR_COMBINED_MULTI=${OPEN_REGIONS_DIR_PAIRS}/timeless_combined_multi
mkdir -p $TIMELESS_DIR_COMBINED_MULTI



CUR_DIR_COMBINED=$TIMELESS_DIR_COMBINED_MULTI/${NUM_CLUSTERS_PROM_ENH}_${NUM_CLUSTERS_PROM_PROM}_${NUM_CLUSTERS_ENH_ENH}
mkdir -p $CUR_DIR_COMBINED_MULTI


qsub -l os=centos7 -l h_vmem=20G -pe smp 1 -cwd -V -m ea \
    -M henriette.miko@mdc-berlin.de -j y \
    -o $CUR_DIR_COMBINED_MULTI/call_matlab_cluster_pairs_combined_out.txt \
    $SCRIPT_DIR/open_regions_subset/open_regions_pairs/call_matlab_cluster_pairs_combined_multi.sh \
    $NUM_MARKS $NUM_FC $SCRIPT_DIR $SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH \
    $SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM $SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH \
    $TIMELESS_DIR_MULTI_PROM_ENH $TIMELESS_DIR_MULTI_PROM_PROM \
    $TIMELESS_DIR_MULTI_ENH_ENH $NUM_CLUSTERS_PROM_ENH $NUM_CLUSTERS_PROM_PROM \
    $NUM_CLUSTERS_ENH_ENH $CUR_DIR_COMBINED_MULTI


exit

