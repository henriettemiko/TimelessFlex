#!/bin/bash

##########
#name:          02b_plot_JAMM_peak_widths.sh
#description:   plot widths of peaks called with JAMM
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


source ../set_variables_hg19.sh

QUALITY_DIR=$OUTPUT_DIR/quality_ChIP

mkdir -p $QUALITY_DIR
cd $QUALITY_DIR

Rscript $SCRIPT_DIR/ChIP/plot_ChIP_JAMM_peak_widths.r $OUTPUT_DIR H3K27ac
Rscript $SCRIPT_DIR/ChIP/plot_ChIP_JAMM_peak_widths.r $OUTPUT_DIR H3K27me3
Rscript $SCRIPT_DIR/ChIP/plot_ChIP_JAMM_peak_widths.r $OUTPUT_DIR H3K4me3
Rscript $SCRIPT_DIR/ChIP/plot_ChIP_JAMM_peak_widths.r $OUTPUT_DIR H3K4me1


exit

