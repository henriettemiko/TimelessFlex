#!/bin/bash

##########
#name:          call_matlab_cluster_pairs_combined_multi.sh
#description:   calls MATLAB for combining models and clustering multi pairs
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 25, 2019
##########


echo $0 started on `hostname` at `date` with parameters $*

NUM_MARKS=$1
NUM_FC=$2
SCRIPT_DIR=$3
SIGNAL_GENERATOR_DIR_PROM_ENH=$4
SIGNAL_GENERATOR_DIR_PROM_PROM=$5
SIGNAL_GENERATOR_DIR_ENH_ENH=$6
MODEL_DIR_PROM_ENH=$7
MODEL_DIR_PROM_PROM=$8
MODEL_DIR_ENH_ENH=$9
CLUSTER_NUM_PROM_ENH=${10}
CLUSTER_NUM_PROM_PROM=${11}
CLUSTER_NUM_ENH_ENH=${12}
CUR_DIR=${13}

matlab -nodisplay -nosplash -nodesktop -r "numTracks=$NUM_MARKS;\
numTimePts=$NUM_FC;signalGeneratorDirPromEnh='$SIGNAL_GENERATOR_DIR_PROM_ENH';\
signalGeneratorDirPromProm='$SIGNAL_GENERATOR_DIR_PROM_PROM';\
signalGeneratorDirEnhEnh='$SIGNAL_GENERATOR_DIR_ENH_ENH';\
modelDirPromEnh='$MODEL_DIR_PROM_ENH';\
modelDirPromProm='$MODEL_DIR_PROM_PROM';\
modelDirEnhEnh='$MODEL_DIR_ENH_ENH';\
numClustersPromEnh=$CLUSTER_NUM_PROM_ENH;\
numClustersPromProm=$CLUSTER_NUM_PROM_PROM;\
numClustersEnhEnh=$CLUSTER_NUM_ENH_ENH;curDir='$CUR_DIR';" < \
    $SCRIPT_DIR/open_regions_subset/open_regions_pairs/\
cluster_pairs_multi_combined.m


exit

