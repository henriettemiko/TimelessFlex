#!/bin/bash

##########
#name:          call_plot_JAMM_peak_widths.sh
#description:   call plot JAMM peak width Rscripts
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          March 28, 2020
##########


OUTPUT_DIR=$1
MARK=$2
SCRIPT_DIR=$3
QUALITY_DIR=$4

echo "#####################"
echo $0 started on `hostname` at `date` with parameters $*
echo "Started processing:"
echo "current mark: $MARK"
echo "#####################"


Rscript $SCRIPT_DIR/ChIP/plot_ChIP_JAMM_peak_widths.r $OUTPUT_DIR $MARK $QUALITY_DIR


exit

