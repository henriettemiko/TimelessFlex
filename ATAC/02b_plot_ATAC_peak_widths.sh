#!/bin/bash

##########
#name:          02b_plot_ATAC_peak_widths.sh
#description:   plot widths of peaks called with MACS2
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


source ../set_variables_hg19.sh

QUALITY_DIR=$OUTPUT_DIR/quality_ATAC

mkdir -p $QUALITY_DIR
cd $QUALITY_DIR


Rscript $SCRIPT_DIR/ATAC/plot_ATAC_MACS2_peak_widths.r $OUTPUT_DIR


exit

