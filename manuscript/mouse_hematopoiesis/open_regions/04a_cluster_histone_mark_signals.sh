#!/bin/bash

##########
#name:          04a_cluster_histone_mark_signals.sh
#description:   calls clustering scripts
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          April 19,2020
##########


source ../set_variables_mm10.sh


NUM_MARKS=4
NUM_FC=5

SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR/Signal_Generator_promoters
SIGNAL_GENERATOR_DIR_ENH=$OPEN_REGIONS_DIR/Signal_Generator_enhancers
TIMELESS_DIR_PROM=$OPEN_REGIONS_DIR/timeless_promoters
mkdir -p $TIMELESS_DIR_PROM
TIMELESS_DIR_ENH=$OPEN_REGIONS_DIR/timeless_enhancers
mkdir -p $TIMELESS_DIR_ENH

START_NUM=2
END_NUM=30

CUR_DIR_ENH=$TIMELESS_DIR_ENH/${START_NUM}_${END_NUM}
mkdir -p $CUR_DIR_ENH
qsub -l os=centos7 -l h_vmem=150G -pe smp 1 -cwd -V -m ea \
    -M henriette.miko@mdc-berlin.de -j y -o $CUR_DIR_ENH/call_matlab_out.txt \
    $SCRIPT_DIR/open_regions/call_matlab.sh $NUM_MARKS $NUM_FC \
    $SCRIPT_DIR $SIGNAL_GENERATOR_DIR_ENH $START_NUM $END_NUM $CUR_DIR_ENH

CUR_DIR_PROM=$TIMELESS_DIR_PROM/${START_NUM}_${END_NUM}
mkdir -p $CUR_DIR_PROM
qsub -l os=centos7 -l h_vmem=150G -pe smp 1 -cwd -V -m ea \
    -M henriette.miko@mdc-berlin.de -j y -o $CUR_DIR_PROM/call_matlab_out.txt \
    $SCRIPT_DIR/open_regions/call_matlab.sh $NUM_MARKS $NUM_FC \
    $SCRIPT_DIR $SIGNAL_GENERATOR_DIR_PROM $START_NUM $END_NUM $CUR_DIR_PROM

exit

