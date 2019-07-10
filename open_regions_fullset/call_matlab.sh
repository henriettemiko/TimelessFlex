#!/bin/bash

##########
#name:          call_matlab.sh
#description:   calls MATLAB for clustering
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 9, 2019
##########


echo $0 started on `hostname` at `date` with parameters $*

NUM_MARKS=$1
NUM_FC=$2
SCRIPT_DIR=$3
SIGNAL_GENERATOR_DIR=$4
START_NUM=$5
END_NUM=$6
CUR_DIR=$7

matlab -nodisplay -nosplash -nodesktop -r "numTracks=$NUM_MARKS;\
    numTimePts=$NUM_FC;endNum=$END_NUM;startNum=$START_NUM;\
    signalGeneratorDir='$SIGNAL_GENERATOR_DIR';curDir='$CUR_DIR';" < \
    $SCRIPT_DIR/open_regions/learn_cluster.m


exit

