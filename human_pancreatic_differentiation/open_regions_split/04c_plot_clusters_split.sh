#!/bin/bash

##########
#name:          04c_plot_clusters_split.sh
#description:   calls plotting of promoter or enhancer clusters for split
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 12, 2019
##########


source ../set_variables_hg19.sh

#promoters


#choose cluster number for promoters here
NUM_CLUSTER_PROM=21


NUM_MARKS=4
NUM_TIME_POINTS=6

SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR_SPLIT/Signal_Generator_promoters
TIMELESS_DIR_PROM=$OPEN_REGIONS_DIR_SPLIT/timeless_promoters


MODEL_DIR_PROM=$TIMELESS_DIR_PROM/2_30
NUM_CLUSTER_DIR_PROM=$TIMELESS_DIR_PROM/2_30/$NUM_CLUSTER_PROM
mkdir -p $NUM_CLUSTER_DIR_PROM
cd $NUM_CLUSTER_DIR_PROM

paste $SIGNAL_GENERATOR_DIR_PROM/allCountsNorm.txt \
    $MODEL_DIR_PROM/classes-${NUM_CLUSTER_PROM}.txt > \
    allCountsNorm_${NUM_CLUSTER_PROM}classes.txt

Rscript $SCRIPT_DIR/open_regions_split/plot_clusters_split.r \
    $NUM_CLUSTER_PROM $NUM_MARKS $NUM_TIME_POINTS $MODEL_DIR_PROM \
    $TIMELESS_DIR_PROM $SIGNAL_GENERATOR_DIR_PROM "promoter"

exit
#enhancers

#choose cluster number for enhancers here
NUM_CLUSTER_ENH=10

SIGNAL_GENERATOR_DIR_ENH=$OPEN_REGIONS_DIR_SPLIT/Signal_Generator_enhancers
TIMELESS_DIR_ENH=$OPEN_REGIONS_DIR_SPLIT/timeless_enhancers

NUM_MARKS=4
NUM_TIME_POINTS=5

MODEL_DIR_ENH=$TIMELESS_DIR_ENH/2_30
NUM_CLUSTER_DIR_ENH=$TIMELESS_DIR_ENH/2_30/$NUM_CLUSTER_ENH
mkdir -p $NUM_CLUSTER_DIR_ENH
cd $NUM_CLUSTER_DIR_ENH

paste $SIGNAL_GENERATOR_DIR_ENH/allCountsNorm.txt \
    $MODEL_DIR_ENH/classes-${NUM_CLUSTER_ENH}.txt > \
    allCountsNorm_${NUM_CLUSTER_ENH}classes.txt

Rscript $SCRIPT_DIR/open_regions_split/plot_clusters_split.r \
    $NUM_CLUSTER_ENH $NUM_MARKS $NUM_TIME_POINTS $MODEL_DIR_ENH \
    $TIMELESS_DIR_ENH $SIGNAL_GENERATOR_DIR_ENH "enhancer"


exit

