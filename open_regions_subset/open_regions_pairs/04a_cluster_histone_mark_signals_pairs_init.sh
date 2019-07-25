#!/bin/bash

##########
#name:          04a_cluster_histone_mark_signals_pairs_init.sh
#description:   calls clustering scripts for pairs init
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 17, 2019
##########


source ../../set_variables_hg19.sh


NUM_MARKS=4
NUM_FC=6 #3 for promoters and 3 for enhancers
START_NUM=2
END_NUM=30

SIGNAL_GENERATOR_DIR_INIT_PROM_ENH=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_init_prom-enh

SIGNAL_GENERATOR_DIR_INIT_PROM_PROM=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_init_prom-prom

SIGNAL_GENERATOR_DIR_INIT_ENH_ENH=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_init_enh-enh


TIMELESS_DIR_INIT_PROM_ENH=$OPEN_REGIONS_DIR_PAIRS/timeless_init_prom-enh
TIMELESS_DIR_INIT_PROM_PROM=$OPEN_REGIONS_DIR_PAIRS/timeless_init_prom-prom
TIMELESS_DIR_INIT_ENH_ENH=$OPEN_REGIONS_DIR_PAIRS/timeless_init_enh-enh

mkdir -p $TIMELESS_DIR_INIT_PROM_ENH
mkdir -p $TIMELESS_DIR_INIT_PROM_PROM
mkdir -p $TIMELESS_DIR_INIT_ENH_ENH


CUR_DIR_INIT_PROM_ENH=$TIMELESS_DIR_INIT_PROM_ENH/${START_NUM}_${END_NUM}
mkdir -p $CUR_DIR_INIT_PROM_ENH

qsub -l os=centos7 -l h_vmem=20G -cwd -V -m ea \
    -M henriette.miko@mdc-berlin.de -j y \
    -o $CUR_DIR_INIT_PROM_ENH/call_matlab_pairs_out.txt \
    $SCRIPT_DIR/open_regions_subset/open_regions_pairs/call_matlab_pairs.sh \
    $NUM_MARKS $NUM_FC $SCRIPT_DIR $SIGNAL_GENERATOR_DIR_INIT_PROM_ENH \
    $START_NUM $END_NUM $CUR_DIR_INIT_PROM_ENH


CUR_DIR_INIT_PROM_PROM=$TIMELESS_DIR_INIT_PROM_PROM/${START_NUM}_${END_NUM}
mkdir -p $CUR_DIR_INIT_PROM_PROM

qsub -l os=centos7 -l h_vmem=20G -cwd -V -m ea \
    -M henriette.miko@mdc-berlin.de -j y \
    -o $CUR_DIR_INIT_PROM_PROM/call_matlab_pairs_out.txt \
    $SCRIPT_DIR/open_regions_subset/open_regions_pairs/call_matlab_pairs.sh \
    $NUM_MARKS $NUM_FC $SCRIPT_DIR $SIGNAL_GENERATOR_DIR_INIT_PROM_PROM \
    $START_NUM $END_NUM $CUR_DIR_INIT_PROM_PROM


CUR_DIR_INIT_ENH_ENH=$TIMELESS_DIR_INIT_ENH_ENH/${START_NUM}_${END_NUM}
mkdir -p $CUR_DIR_INIT_ENH_ENH

qsub -l os=centos7 -l h_vmem=20G -cwd -V -m ea \
    -M henriette.miko@mdc-berlin.de -j y \
    -o $CUR_DIR_INIT_ENH_ENH/call_matlab_pairs_out.txt \
    $SCRIPT_DIR/open_regions_subset/open_regions_pairs/call_matlab_pairs.sh \
    $NUM_MARKS $NUM_FC $SCRIPT_DIR $SIGNAL_GENERATOR_DIR_INIT_ENH_ENH \
    $START_NUM $END_NUM $CUR_DIR_INIT_ENH_ENH


exit

