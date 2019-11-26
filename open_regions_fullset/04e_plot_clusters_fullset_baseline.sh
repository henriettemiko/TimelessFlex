#!/bin/bash

##########
#name:          04e_plot_clusters_fullset_baseline.sh
#description:   calls plotting baseline approaches for fullset
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          November 21, 2019
##########


source ../set_variables_hg19.sh

#promoters


#choose cluster number for promoters here
NUM_CLUSTER_PROM=14


NUM_MARKS=4
NUM_TIME_POINTS=5

SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR_FULL/Signal_Generator_promoters
TIMELESS_DIR_PROM=$OPEN_REGIONS_DIR_FULL/timeless_promoters


MODEL_DIR_PROM=$TIMELESS_DIR_PROM/2_30
NUM_CLUSTER_DIR_PROM=$TIMELESS_DIR_PROM/2_30/$NUM_CLUSTER_PROM
mkdir -p $NUM_CLUSTER_DIR_PROM
cd $NUM_CLUSTER_DIR_PROM

#paste $SIGNAL_GENERATOR_DIR_PROM/allCountsNorm.txt \
#    $MODEL_DIR_PROM/classes-${NUM_CLUSTER_PROM}.txt > \
#    allCountsNorm_${NUM_CLUSTER_PROM}classes.txt


Rscript $SCRIPT_DIR/open_regions_fullset/baseline_kmeans_Zscore.r \
    $NUM_CLUSTER_PROM $NUM_MARKS $NUM_TIME_POINTS $MODEL_DIR_PROM \
    $TIMELESS_DIR_PROM $SIGNAL_GENERATOR_DIR_PROM "promoter" "prom" 


Rscript $SCRIPT_DIR/open_regions_fullset/baseline_kmeans_normcounts.r \
    $NUM_CLUSTER_PROM $NUM_MARKS $NUM_TIME_POINTS $MODEL_DIR_PROM \
    $TIMELESS_DIR_PROM $SIGNAL_GENERATOR_DIR_PROM "promoter" "prom" 

Rscript $SCRIPT_DIR/open_regions_fullset/baseline_kmeans_normfold.r \
    $NUM_CLUSTER_PROM $NUM_MARKS $NUM_TIME_POINTS $MODEL_DIR_PROM \
    $TIMELESS_DIR_PROM $SIGNAL_GENERATOR_DIR_PROM "promoter" "prom"
 
#exit
#enhancers

#choose cluster number for enhancers here
NUM_CLUSTER_ENH=15

SIGNAL_GENERATOR_DIR_ENH=$OPEN_REGIONS_DIR_FULL/Signal_Generator_enhancers
TIMELESS_DIR_ENH=$OPEN_REGIONS_DIR_FULL/timeless_enhancers

NUM_MARKS=4
NUM_TIME_POINTS=5

MODEL_DIR_ENH=$TIMELESS_DIR_ENH/2_30
NUM_CLUSTER_DIR_ENH=$TIMELESS_DIR_ENH/2_30/$NUM_CLUSTER_ENH
mkdir -p $NUM_CLUSTER_DIR_ENH
cd $NUM_CLUSTER_DIR_ENH

#paste $SIGNAL_GENERATOR_DIR_ENH/allCountsNorm.txt \
#    $MODEL_DIR_ENH/classes-${NUM_CLUSTER_ENH}.txt > \
#    allCountsNorm_${NUM_CLUSTER_ENH}classes.txt

Rscript $SCRIPT_DIR/open_regions_fullset/baseline_kmeans_Zscore.r \
    $NUM_CLUSTER_ENH $NUM_MARKS $NUM_TIME_POINTS $MODEL_DIR_ENH \
    $TIMELESS_DIR_ENH $SIGNAL_GENERATOR_DIR_ENH "enhancer" "enh"

Rscript $SCRIPT_DIR/open_regions_fullset/baseline_kmeans_normcounts.r \
    $NUM_CLUSTER_ENH $NUM_MARKS $NUM_TIME_POINTS $MODEL_DIR_ENH \
    $TIMELESS_DIR_ENH $SIGNAL_GENERATOR_DIR_ENH "enhancer" "enh"

Rscript $SCRIPT_DIR/open_regions_fullset/baseline_kmeans_normfold.r \
    $NUM_CLUSTER_ENH $NUM_MARKS $NUM_TIME_POINTS $MODEL_DIR_ENH \
    $TIMELESS_DIR_ENH $SIGNAL_GENERATOR_DIR_ENH "enhancer" "enh"

exit

