#!/bin/bash

##########
#name:          call_matlab_cluster_multi_pairs_combined.sh
#description:   calls MATLAB for combining models and clustering all pairs
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 25, 2019
##########


echo $0 started on `hostname` at `date` with parameters $*

NUM_MARKS=$1
NUM_FC=$2
SCRIPT_DIR=$3
SIGNAL_GENERATOR_DIR_INIT_PROM_ENH=$4
SIGNAL_GENERATOR_DIR_INIT_PROM_PROM=$5
SIGNAL_GENERATOR_DIR_INIT_ENH_ENH=$6
SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH=$7
SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM=$8
SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH=$9
MODEL_DIR_PROM_ENH=${10}
MODEL_DIR_PROM_PROM=${11}
MODEL_DIR_ENH_ENH=${12}
CLUSTER_NUM_PROM_ENH=${13}
CLUSTER_NUM_PROM_PROM=${14}
CLUSTER_NUM_ENH_ENH=${15}
CUR_DIR=${16}

matlab -nodisplay -nosplash -nodesktop -r "numTracks=$NUM_MARKS;\
numTimePts=$NUM_FC;\
signalGeneratorDirPromEnhInit='$SIGNAL_GENERATOR_DIR_INIT_PROM_ENH';\
signalGeneratorDirPromPromInit='$SIGNAL_GENERATOR_DIR_INIT_PROM_PROM';\
signalGeneratorDirEnhEnhInit='$SIGNAL_GENERATOR_DIR_INIT_ENH_ENH';\
signalGeneratorDirPromEnhMulti='$SIGNAL_GENERATOR_DIR_MULTI_PROM_ENH';\
signalGeneratorDirPromPromMulti='$SIGNAL_GENERATOR_DIR_MULTI_PROM_PROM';\
signalGeneratorDirEnhEnhMulti='$SIGNAL_GENERATOR_DIR_MULTI_ENH_ENH';\
modelDirPromEnh='$MODEL_DIR_PROM_ENH';\
modelDirPromProm='$MODEL_DIR_PROM_PROM';\
modelDirEnhEnh='$MODEL_DIR_ENH_ENH';\
numClustersPromEnh=$CLUSTER_NUM_PROM_ENH;\
numClustersPromProm=$CLUSTER_NUM_PROM_PROM;\
numClustersEnhEnh=$CLUSTER_NUM_ENH_ENH;curDir='$CUR_DIR';" < \
    $SCRIPT_DIR/open_regions_subset/open_regions_pairs/cluster_pairs_initmulti_combined.m

exit

