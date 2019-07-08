#!/bin/bash

source ../set_variables_hg19.sh


QUALITY_DIR=$OUTPUT_DIR/quality_ChIP

mkdir -p $QUALITY_DIR
cd $QUALITY_DIR


#JAMM peaks
Rscript $SCRIPT_DIR/ChIP/plot_ChIP_JAMM_peak_widths.r $OUTPUT_DIR H3K27ac
Rscript $SCRIPT_DIR/ChIP/plot_ChIP_JAMM_peak_widths.r $OUTPUT_DIR H3K27me3
Rscript $SCRIPT_DIR/ChIP/plot_ChIP_JAMM_peak_widths.r $OUTPUT_DIR H3K4me3
Rscript $SCRIPT_DIR/ChIP/plot_ChIP_JAMM_peak_widths.r $OUTPUT_DIR H3K4me1


#MACS2 peaks
Rscript $SCRIPT_DIR/ChIP/plot_ChIP_MACS2_peak_widths.r $OUTPUT_DIR H3K27ac
Rscript $SCRIPT_DIR/ChIP/plot_ChIP_MACS2_peak_widths.r $OUTPUT_DIR H3K27me3
Rscript $SCRIPT_DIR/ChIP/plot_ChIP_MACS2_peak_widths.r $OUTPUT_DIR H3K4me3
Rscript $SCRIPT_DIR/ChIP/plot_ChIP_MACS2_peak_widths.r $OUTPUT_DIR H3K4me1


exit
