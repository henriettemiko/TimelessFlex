#!/bin/bash

##########
#name:          04c_plot_clusters_pairs_init.sh
#description:   calls plotting of clusters for pairs init
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 17, 2019
##########


source ../../set_variables_hg19.sh


#choose cluster number for init prom-enh pairs here
NUM_CLUSTER_INIT_PROM_ENH=10

#choose cluster number for init prom-prom pairs here
NUM_CLUSTER_INIT_PROM_PROM=4

#choose cluster number for init enh-enh pairs here
NUM_CLUSTER_INIT_ENH_ENH=12


NUM_MARKS=4
NUM_TIME_POINTS=4


#init prom-enh

SIGNAL_GENERATOR_DIR_INIT_PROM_ENH=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_init_prom-enh
TIMELESS_DIR_INIT_PROM_ENH=$OPEN_REGIONS_DIR_PAIRS/timeless_init_prom-enh

MODEL_DIR_INIT_PROM_ENH=$TIMELESS_DIR_INIT_PROM_ENH/2_30
NUM_CLUSTER_DIR_INIT_PROM_ENH=$TIMELESS_DIR_INIT_PROM_ENH/2_30/\
$NUM_CLUSTER_INIT_PROM_ENH
mkdir -p $NUM_CLUSTER_DIR_INIT_PROM_ENH

cd $NUM_CLUSTER_DIR_INIT_PROM_ENH

paste $SIGNAL_GENERATOR_DIR_INIT_PROM_ENH/allCountsNorm.txt \
    $MODEL_DIR_INIT_PROM_ENH/\
classes-${NUM_CLUSTER_INIT_PROM_ENH}.txt > \
    allCountsNorm_${NUM_CLUSTER_INIT_PROM_ENH}classes.txt

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
plot_clusters_pairs.r \
    $NUM_CLUSTER_INIT_PROM_ENH $NUM_MARKS $NUM_TIME_POINTS \
    $MODEL_DIR_INIT_PROM_ENH $TIMELESS_DIR_INIT_PROM_ENH \
    $SIGNAL_GENERATOR_DIR_INIT_PROM_ENH "init_prom-enh"


#init prom-prom

SIGNAL_GENERATOR_DIR_INIT_PROM_PROM=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_init_prom-prom
TIMELESS_DIR_INIT_PROM_PROM=$OPEN_REGIONS_DIR_PAIRS/timeless_init_prom-prom

MODEL_DIR_INIT_PROM_PROM=$TIMELESS_DIR_INIT_PROM_PROM/2_30
NUM_CLUSTER_DIR_INIT_PROM_PROM=$TIMELESS_DIR_INIT_PROM_PROM/2_30/\
$NUM_CLUSTER_INIT_PROM_PROM
mkdir -p $NUM_CLUSTER_DIR_INIT_PROM_PROM

cd $NUM_CLUSTER_DIR_INIT_PROM_PROM

paste $SIGNAL_GENERATOR_DIR_INIT_PROM_PROM/allCountsNorm.txt \
    $MODEL_DIR_INIT_PROM_PROM/\
classes-${NUM_CLUSTER_INIT_PROM_PROM}.txt > \
    allCountsNorm_${NUM_CLUSTER_INIT_PROM_PROM}classes.txt

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
plot_clusters_pairs.r \
    $NUM_CLUSTER_INIT_PROM_PROM $NUM_MARKS $NUM_TIME_POINTS \
    $MODEL_DIR_INIT_PROM_PROM $TIMELESS_DIR_INIT_PROM_PROM \
    $SIGNAL_GENERATOR_DIR_INIT_PROM_PROM "init_prom-prom"


#init enh-enh

SIGNAL_GENERATOR_DIR_INIT_ENH_ENH=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_init_enh-enh
TIMELESS_DIR_INIT_ENH_ENH=$OPEN_REGIONS_DIR_PAIRS/timeless_init_enh-enh

MODEL_DIR_INIT_ENH_ENH=$TIMELESS_DIR_INIT_ENH_ENH/2_30
NUM_CLUSTER_DIR_INIT_ENH_ENH=$TIMELESS_DIR_INIT_ENH_ENH/2_30/\
$NUM_CLUSTER_INIT_ENH_ENH
mkdir -p $NUM_CLUSTER_DIR_INIT_ENH_ENH

cd $NUM_CLUSTER_DIR_INIT_ENH_ENH

paste $SIGNAL_GENERATOR_DIR_INIT_ENH_ENH/allCountsNorm.txt \
    $MODEL_DIR_INIT_ENH_ENH/\
classes-${NUM_CLUSTER_INIT_ENH_ENH}.txt > \
    allCountsNorm_${NUM_CLUSTER_INIT_ENH_ENH}classes.txt

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
plot_clusters_pairs.r \
    $NUM_CLUSTER_INIT_ENH_ENH $NUM_MARKS $NUM_TIME_POINTS \
    $MODEL_DIR_INIT_ENH_ENH $TIMELESS_DIR_INIT_ENH_ENH \
    $SIGNAL_GENERATOR_DIR_INIT_ENH_ENH "init_enh-enh"


exit

