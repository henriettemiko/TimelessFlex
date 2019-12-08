#!/bin/bash

##########
#name:          call_plot_ATAC_MACS2_peak_widths_D8.sh
#description:   call plot ATAC MACS2 peak width Rscripts
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 27, 2019
##########


OUTPUT_DIR=$1
SCRIPT_DIR=$2
QUALITY_DIR=$3

echo "#####################"
echo $0 started on `hostname` at `date` with parameters $*
echo "Started processing:"
echo "#####################"


Rscript $SCRIPT_DIR/ATAC/plot_ATAC_MACS2_peak_widths_D8.r $OUTPUT_DIR $QUALITY_DIR


exit

