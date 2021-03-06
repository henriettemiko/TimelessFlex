#!/bin/bash

##########
#name:          call_matlab_multi_pairs.sh
#description:   calls MATLAB for clustering pairs multi
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 17, 2019
##########


echo $0 started on `hostname` at `date` with parameters $*

NUM_MARKS=$1
NUM_FC=$2
SCRIPT_DIR=$3
SIGNAL_GENERATOR_DIR=$4
MODEL_DIR=$5
CLUSTER_NUM=$6
CUR_DIR=$7

matlab -nodisplay -nosplash -nodesktop -r "numTracks=$NUM_MARKS;\
numTimePts=$NUM_FC;signalGeneratorDir='$SIGNAL_GENERATOR_DIR';\
numClusters=$CLUSTER_NUM;modelDir='$MODEL_DIR';curDir='$CUR_DIR';" < \
    $SCRIPT_DIR/open_regions_subset/open_regions_pairs/cluster_multi_pairs.m


exit

