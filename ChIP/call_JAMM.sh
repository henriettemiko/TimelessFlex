#!/bin/bash

##########
#name:          call_JAMM.sh
#description:   call peaks with JAMM for each time point
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


JAMM_DIR=$1
MARK=$2
TIME=$3
OUTPUT_DIR=$4
CHR_SIZES=$5

echo "#####################"
echo $0 started on `hostname` at `date` with parameters $*
echo "Started processing:"
echo "current histone marks: $MARK"
echo "current time point: $TIME"
echo "#####################"


mkdir -p $JAMM_DIR/tmp

bash JAMM.sh -s $OUTPUT_DIR/ChIP/$MARK/$TIME/mapping \
    -c $OUTPUT_DIR/ChIP/Input/$TIME/mapping -o $JAMM_DIR -g $CHR_SIZES \
    -r window -e auto -b 150 -T $JAMM_DIR/tmp


exit

